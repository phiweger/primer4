import streamlit as st


def warn(method, params, amplicon_len_min, amplicon_len_max):
    lu = {'sanger': 'size_range_PCR', 'qpcr': 'size_range_qPCR', 'mrna': 'size_range_mRNA'}
    # Overwrite params from config
    original = params[lu[method]]
    params[lu[method]] = [amplicon_len_min, amplicon_len_max]
    # import warnings

    if not all(i==j for i, j in zip(original, [amplicon_len_min, amplicon_len_max])):
        # https://emojifinder.com/fear
        st.warning('Amplicon values outside of typical range')
        # icon='😱'
    return None