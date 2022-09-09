from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st


def prepare_data_for_vis():
    # Data for plotting
    # https://docs.streamlit.io/library/api-reference/charts/st.pyplot
    return np.random.normal(1, 1, size=300)


class Beauty():
    def __init__(self, data):
        self.data = data
        
    def plot(self, width=10, height=4, dpi=300):
        fig, ax = plt.subplots(figsize=(width, height))
        # fig.set_size_inches(1, 1)
        fig.set_dpi(dpi)
        ax.hist(self.data, bins=30)

        # Oh FFS, there's always something:
        # discuss.streamlit.io/t/cannot-change-matplotlib-figure-size/10295/8
        # TL;DR: st.pyplot(fig) would be ideal but does not scale img correctly
        
        # Workaround for streamlit to respect figure size in matplotlib objects
        buf = BytesIO()
        fig.savefig(buf, format="png")
        # Workaround for streamlit to center the image
        col1, col2, col3 = st.columns([1, 1000, 1])
        with col1:
            st.write('')
        with col2:
            st.image(buf)
        with col3:
            st.write('')    

        return None
