"""
modules/preparer.py
"""

import os
import glob

import seekr.modules.structures as structures
import seekr.modules.engines.structures as engines_structures

def cleanse_outputs_stage(
        model: "structures.Seekr_model",
        anchor_or_unpartitioned: engines_structures.MD_unpartitioned_settings \
            | structures.Anchor, 
        stage_index: int
        ) -> None:
    """
    Delete all simulation outputs from unpartitioned settings or anchor. If, say, a 
    reporter generates files, they can be deleted here.
    """
    stage = model.stages[stage_index]
    if stage.type == "md_stage":
        #full_directory = os.path.join(model.directory, anchor_or_unpartitioned.directory, 
        #                                stage.get_directory_name())
        #if os.path.exists(full_directory):
        #    file_glob = os.path.join(full_directory, "*.WHATEVER")
        #    for myfile in glob.glob(file_glob):
        #        os.remove(myfile)
        pass
