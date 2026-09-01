from .widgets import PageDescriptor

# Compatibility module: chart projections live on the Charts page; PyMOL
# controls have their own page in the integrated workflow.
PAGE = PageDescriptor("Charts", "Inspect authoritative chart projections and their explicit units.")

__all__ = ["PAGE"]
