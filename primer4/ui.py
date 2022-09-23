'''
hier ist das vrex aggrid, du brauchst def set_grid_options_variants_for_vrex(self, df) und def draw_grid_variants_for_vrex(self, df)
'''



'''

        # Examples:
        # https://pablocfonseca-streamlit-aggrid-examples-example-jyosi3.streamlitapp.com/
        from st_aggrid import AgGrid, GridOptionsBuilder
        from st_aggrid.shared import GridUpdateMode

        # https://streamlit-aggrid.readthedocs.io/en/docs/AgGrid.html
        # 
        #response = AgGrid(df, editable=True, fit_columns_on_grid_load=True, update_mode=GridUpdateMode.NO_UPDATE, reload_data=False)
        # st.write(response['data'])

        st.write("double click fields to modify and select for import")

        options = GridOptionsBuilder.from_dataframe(
                        df,
                        enableRowGroup=False,
                        enableValue=False,
                        enablePivot=False
                        )
        options.configure_column('Amplicon', 
                            editable=True,
                            cellEditor="agRichSelectCellEditor",
                            cellEditorPopup=False
                            )
        variant_upload_selection = AgGrid(
                                df,
                                enable_enterprise_modules=True,
                                gridOptions=options.build(),
                                theme="balham",
                                update_mode=GridUpdateMode.MODEL_CHANGED,
                                allow_unsafe_jscode=True,
                                height=300)
        pd.json_normalize(variant_upload_selection["selected_rows"])
'''




################################################################
# general packages
import pandas as pd
import time

################################################################
# streamlit packages
import streamlit as st

import streamlit_modal as modal
# TODO: streamlit will release native modals soon

from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import \
    GridUpdateMode, \
    JsCode

################################################################
# vrex packages

from Output.GEPADO import GEPADO_variant_writer

from Streamlit_components import \
    streamlit_labels as sl

from Streamlit_components import \
    session_names

from Settings import \
    variant_attributes as va

################################################################


class import_grid:
    """
    This class takes a dataframe and displays a grid for
    selection and modification of variant attributes.
    Mandatory fields are checked.
    Selected fields are updated with hgnc id.
    Draws an import button and imports variants to vrex. 
    
    parameters:
        df_snvs_for_vrex: dataframe with variants, ready for import

    """ 
    def __init__(self, df_variants_for_vrex, variant_type, current_analysis_id):
        self.df_variants_for_vrex = self.clean_up_df(df_variants_for_vrex)
        self.current_analysis_id = current_analysis_id
        self.variant_type = variant_type
        
    def clean_up_df(self, df):
        df.fillna("",inplace=True)
        self.na_strings = ["nan", "NA", "", "None"]
        for item in self.na_strings:
            df = df.replace(item, "", regex=True)
        return df

    def set_grid_options_variants_for_vrex(self, df):

        self.options = GridOptionsBuilder.from_dataframe(
                        df,
                        enableRowGroup=False,
                        enableValue=False,
                        enablePivot=False
                        )
        self.options.configure_side_bar()
        self.options.configure_selection("multiple", use_checkbox=True)
    
        if self.variant_type == "snv":
            # configure freetext
            for attribute in va.snv_free_text_attributes:
                self.options.configure_column(attribute, editable=True)

            # configure dropdowns
            for key in va.snv_dropdown_attributes:
                self.options.configure_column(key, 
                    editable=True,
                    cellEditor="agRichSelectCellEditor",
                    cellEditorParams={"values":va.snv_dropdown_attributes[key]},
                    cellEditorPopup=False
                    )

            # color coding mandatory cells
            for attribute in va.snv_mandatory_attributes:
                self.options.configure_column(attribute, cellStyle=JsCode(va.cellstyle_jscode))

        if self.variant_type == "cnv":
            # configure freetext
            for attribute in va.cnv_free_text_attributes:
                self.options.configure_column(attribute, editable=True)
            
            # configure dropdowns
            for key in va.cnv_dropdown_attributes:
                self.options.configure_column(key, 
                    editable=True,
                    cellEditor="agRichSelectCellEditor",
                    cellEditorParams={"values":va.cnv_dropdown_attributes[key]},
                    cellEditorPopup=False
                    )

            # color coding mandatory cells
            for attribute in va.cnv_mandatory_attributes:
                self.options.configure_column(attribute, cellStyle=JsCode(va.cellstyle_jscode))

        # "select all"-checkbox on first column head
        self.options.configure_column("Befundnummer", headerCheckboxSelection = True)

    def draw_grid_variants_for_vrex(self, df):
        st.write("double click fields to modify and select for import")

        variant_upload_selection = AgGrid(
                                df,
                                enable_enterprise_modules=True,
                                gridOptions=self.options.build(),
                                theme="balham",
                                update_mode=GridUpdateMode.MODEL_CHANGED,
                                allow_unsafe_jscode=True,
                                height=300)
        self.df_selected_variants = pd.json_normalize(variant_upload_selection["selected_rows"])
        self.df_selected_variants = self.update_hgnc_id(self.df_selected_variants)
        
    
    def update_hgnc_id(self, df):
        # update hgnc id
        if self.variant_type == "snv":
            for idx, row in df.iterrows():
                if row["snv_hgnc_symbol"] in self.na_strings:
                    df.at[idx, "snv_hgnc_id"] = ""
                else:
                    df.at[idx, "snv_hgnc_id"] = va.hgnc_dict[row["snv_hgnc_symbol"]].replace("HGNC:", "")
        return df
    
    def check_mandatory_fields(self):
        if self.variant_type == "snv":
            mandatory_attributes = va.snv_mandatory_attributes
        elif self.variant_type == "cnv":
            mandatory_attributes = va.cnv_mandatory_attributes
        
        self.mandatory_fields_ok = True
        for attribute in mandatory_attributes:
            for idx, row in self.df_selected_variants.iterrows():
                if row[attribute] in self.na_strings:
                    self.mandatory_fields_ok = False
        if not self.mandatory_fields_ok:
            st.warning(sl.manual_import_mandatory_fields_warning)    

    def add_row(self, df_variants):

        df_append = pd.DataFrame(columns=df_variants.columns)
        df_append.at[0, "Befundnummer"] = self.current_analysis_id

        df_variants = df_variants.append(df_append, ignore_index = True)
        return df_variants

    def import_variant_collection(self, db_connector):
        self.df_selected_variants = st.session_state[session_names.df_manual_variant_collection]

        if st.button(sl.manual_import_import_variants, key=self.variant_type+"_import"):
            if self.df_selected_variants.empty:
                st.warning(sl.manual_import_no_variants)
            else:
                gepado_variant_writer = GEPADO_variant_writer.variant_writer(
                    db_connector,
                    self.current_analysis_id
                )
                # prepare sql for import
                sql_snippets = gepado_variant_writer.create_sql_command_for_export(self.df_selected_variants, self.variant_type)
                
                # some basic logging
                for idx, row in self.df_selected_variants.iterrows():
                    print("importing:\n", row)



                with st.spinner(text="importing..."):
                    for snippet in sql_snippets:
                        gepado_variant_writer.execute_vrex_sql(snippet)
                    time.sleep(2)
                    st.write("import sucessful!")
                    
                    st.experimental_rerun()