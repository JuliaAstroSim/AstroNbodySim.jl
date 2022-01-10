"""
$(TYPEDEF)

Controls accuracy of tree method

$(TYPEDFIELDS)
"""
struct TreeSimConfig
    "If the view angle of tree node that relative to the sink point is too large, open this node. In radian unit. Default is `0.1`"
    TreeOpenAngle::Float64
    "While opening the tree node, take the last acceleration into account. Improves accuracy at high redshift cosmology. Default is `0.025`"
    ErrTolAcc::Float64
end

function TreeSimConfig(;
        TreeOpenAngle = 0.1, # rad
        ErrTolAcc = 0.025,
    )
    return TreeSimConfig(
        TreeOpenAngle,
        ErrTolAcc,
    )
end

struct OctreeData{A, U, Len, Len_1, I, F, POS, VEL, MASS, B, Ext}
    treesimconfig::TreeSimConfig
    tree::Octree{A, U, Len, Len_1, I, F, POS, VEL, MASS, B, Ext}
end