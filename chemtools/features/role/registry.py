from __future__ import annotations

REGISTRY = {
    "global": [
        {"name": "global.MW", "type": "float", "default": 0.0},
        {"name": "global.logP", "type": "float", "default": 0.0},
        {"name": "global.TPSA", "type": "float", "default": 0.0},
        {"name": "global.rotB", "type": "int", "default": 0},
        {"name": "global.HBD", "type": "int", "default": 0},
        {"name": "global.HBA", "type": "int", "default": 0},
        {"name": "global.aromatic_rings", "type": "int", "default": 0},
        {"name": "global.heteroatoms", "type": "int", "default": 0},
    ],
    "amine": [
        {"name": "amine.present", "type": "int", "default": 0},
        {"name": "amine.class_ps3", "type": "int", "default": -1},
        {"name": "amine.aniline_flag", "type": "int", "default": 0},
        {"name": "amine.alpha_branch", "type": "int", "default": 0},
        {"name": "amine.formal_charge", "type": "int", "default": 0},
        {"name": "amine.is_aromatic_N", "type": "int", "default": 0},
        {"name": "amine.aryl_neighbors", "type": "int", "default": 0},
        {"name": "amine.alpha_heavy_neighbors", "type": "int", "default": 0},
        {"name": "amine.h_count_on_N", "type": "int", "default": 0},
    ],
    "alcohol": [
        {"name": "alcohol.present", "type": "int", "default": 0},
        {"name": "alcohol.class_ps3", "type": "int", "default": -1},
        {"name": "alcohol.benzylic", "type": "int", "default": 0},
        {"name": "alcohol.allylic", "type": "int", "default": 0},
        {"name": "alcohol.phenol_flag", "type": "int", "default": 0},
        {"name": "alcohol.alpha_degree", "type": "int", "default": -1},
        {"name": "alcohol.alpha_heavy_neighbors", "type": "int", "default": 0},
        {"name": "alcohol.aryl_neighbors", "type": "int", "default": 0},
    ],
    "aryl_halide": [
        {"name": "aryl_halide.present", "type": "int", "default": 0},
        {"name": "aryl_halide.halide", "type": "int", "default": 0},  # 1:F 2:Cl 3:Br 4:I 0:UNK
        {"name": "aryl_halide.ortho_block", "type": "int", "default": 0},
        {"name": "aryl_halide.ipso_degree", "type": "int", "default": 0},
        {"name": "aryl_halide.para_EWG", "type": "int", "default": 0},
        {"name": "aryl_halide.heteroaryl", "type": "int", "default": 0},
        {"name": "aryl_halide.ring_aromatic_count", "type": "int", "default": 0},
        {"name": "aryl_halide.triflate_flag", "type": "int", "default": 0},
    ],
    "fingerprints": {
        "amine": {"bits": 512, "radius": 2},
        "alcohol": {"bits": 512, "radius": 2},
        "aryl_halide": {"bits": 512, "radius": 2},
    },
}
