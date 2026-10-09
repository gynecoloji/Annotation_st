"""annotation_st — config-driven cell annotation for Xenium spatial transcriptomics.

Stages (each also a library function):

* ``qc``       — :func:`annotation_st.stages.qc.qc_filter`
* ``lineage``  — :func:`annotation_st.stages.lineage.annotate_lineage`
* ``caf``      — :func:`annotation_st.stages.caf.annotate_caf`
* ``csc``      — :func:`annotation_st.stages.csc.call_csc`
* ``immune``   — :func:`annotation_st.stages.immune.run_immune`
* ``merge``    — :func:`annotation_st.stages.merge.merge_labels`
* ``spatial``  — :func:`annotation_st.stages.spatial.build_spatial_graph` and friends
"""
from .config import PipelineConfig, default_config, dump_config, load_config, preset_config

# Bumped automatically by release-please on each release — do not edit by hand.
__version__ = "0.2.0"  # x-release-please-version
__all__ = ["PipelineConfig", "default_config", "preset_config", "dump_config", "load_config", "__version__"]
