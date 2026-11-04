"""Reaction on Streamlit Community Cloud: ``streamlit run streamlit_app.py``.

Two windows, as in the desktop app: Input (the scenario card, this address) and Results (``/results``),
sharing one session through the private ``?w=`` id in the address; and the Sources page.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import streamlit as st  # noqa: E402

from reaction.web.app import input_page, results_page, sources_page  # noqa: E402

page = st.navigation([st.Page(input_page, title="Input", url_path="input", default=True),
                      st.Page(results_page, title="Results", url_path="results"),
                      st.Page(sources_page, title="Sources", url_path="sources")], position="hidden")
page.run()
