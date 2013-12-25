

/*
# Construct from Gmsh .msh formatted input stream.
function TriMesh(is::IO, subdiv_iters::Integer,
                 integration_rel_err::R = 1e-12, integration_abs_err::R = 1e-12,
                 load_phys_reg_tags::Bool = false, load_geom_ent_tags::Bool = false)

  const base_pts_by_nodenum = read_gmsh_nodes(is)

  # Read to the beginning of the elements section.
  if !read_through_line_starting(is, "\$Elements")
    error("Could not find beginning of elements section in mesh file.")
  end

  # Estimate number of input mesh triangles from the number of elements declared and the number
  # of skipped lower order elements (points and lines) before the first triangle element.
  const decl_num_els = uint64(readline(is)) # This count can include unwanted lower order elements.
  const line, lines_read = read_until(is, is_polytope_or_endmarker)
  if line == "" || beginswith(line, "\$EndElements")
    error("No elements found in \$Elements section.")
  end
  const est_base_tris = decl_num_els - (lines_read - 1)

  const base_tris_iter = GmshElementsIter(line, is)

  TriMesh(base_pts_by_nodenum,
          base_tris_iter,
          est_base_tris,
          subdiv_iters,
          integration_rel_err, integration_abs_err,
          load_phys_reg_tags, load_geom_ent_tags)
end





# ----------------------------------
# Gmsh msh format reading support

typealias ElTypeNum Uint8
eltypenum(i::Integer) = if i>=1 && i<=typemax(ElTypeNum) convert(ElTypeNum, i) else error("invalid element type: $i") end
eltypenum(s::String) = eltypenum(int(s))

# Gmsh triangle type codes
const ELTYPE_3_NODE_TRIANGLE = eltypenum(2)
# lower order elements, which can be ignored
const ELTYPE_POINT = eltypenum(15)
const ELTYPE_2_NODE_LINE = eltypenum(1)
const ELTYPE_3_NODE_LINE = eltypenum(8)
const ELTYPE_4_NODE_LINE = eltypenum(26)
const ELTYPE_5_NODE_LINE = eltypenum(27)
const ELTYPE_6_NODE_LINE = eltypenum(28)

function is_3_node_triangle_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_3_NODE_TRIANGLE
end

function is_lower_order_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_POINT ||
  el_type == ELTYPE_2_NODE_LINE ||
  el_type == ELTYPE_3_NODE_LINE ||
  el_type == ELTYPE_4_NODE_LINE ||
  el_type == ELTYPE_5_NODE_LINE ||
  el_type == ELTYPE_6_NODE_LINE
end

const TOKN_NODELINE_ELNUM = 1
const TOKN_NODELINE_POINT1 = 2
const TOKN_NODELINE_POINT2 = 3
const TOKN_NODELINE_POINT3 = 4

# Gmsh element line format:
# elm-number elm-type number-of-tags < tag > ... vert-number-list
const TOKN_ELLINE_ELTYPE = 2
const TOKN_ELLINE_NUMTAGS = 3


# Gmsh base elements iterator

immutable GmshElementsIter
  first_base_tri_line::String
  is::IO # stream should be positioned after first base element line
end

import Base.start, Base.done, Base.next
function start(gmsh_els_iter::GmshElementsIter)
  base_tri_from_gmsh_el(split(gmsh_els_iter.first_base_tri_line, ' '))
end

function done(gmsh_els_iter::GmshElementsIter, next_base_tri)
  next_base_tri == nothing
end

function next(gmsh_els_iter::GmshElementsIter, next_base_tri)
  while true
    const line = readline(gmsh_els_iter.is)
    if beginswith(line, "\$EndElements")
      return (next_base_tri, nothing)
    elseif line == ""
      error("No EndElements marker found.")
    else
      const toks = split(line, ' ')
      const el_type = eltypenum(toks[TOKN_ELLINE_ELTYPE])
      if is_lower_order_el_type(el_type)
        continue
      elseif is_3_node_triangle_el_type(el_type)
        return (next_base_tri, base_tri_from_gmsh_el(toks))
      else
        error("Element type $el_type is not supported.")
      end
    end
  end
end

const EMPTY_OTHER_TAGS = Array(Tag,0)

function base_tri_from_gmsh_el(toks::Array{String,1})
  const num_tags = int(toks[TOKN_ELLINE_NUMTAGS])
  assert(num_tags >= 2)
  const physreg_tag = tag(int64(toks[TOKN_ELLINE_NUMTAGS+1]))
  const geoment_tag = tag(int64(toks[TOKN_ELLINE_NUMTAGS+2]))
  const other_tags = num_tags >= 3 ? map(t -> tag(int64(t)), toks[TOKN_ELLINE_NUMTAGS+3:TOKN_ELLINE_NUMTAGS+num_tags]) :
                                     EMPTY_OTHER_TAGS
  const point_nums = let last_tag_tokn = TOKN_ELLINE_NUMTAGS + num_tags;
    (nodenum(toks[last_tag_tokn + 1]),
     nodenum(toks[last_tag_tokn + 2]),
     nodenum(toks[last_tag_tokn + 3]))
  end

  BaseTri(point_nums,
          physreg_tag,
          geoment_tag,
          other_tags)
end


###############################################
# Mesh File Reading Utilities

function read_gmsh_nodes(is::IO)
  if !read_through_line_starting(is, "\$Nodes")
    error("Nodes section not found in mesh input file.")
  else
    # Next line should be the vert count.
    const count = int(readline(is))
    const pts = Array(Point, count)
    l = readline(is)
    while !beginswith(l, "\$EndNodes") && l != ""
      const toks = split(l, ' ')
      const pt = Point(convert(R, float64(toks[TOKN_NODELINE_POINT1])), convert(R, float64(toks[TOKN_NODELINE_POINT2])))
      if length(toks) >= TOKN_NODELINE_POINT3 && float64(toks[TOKN_NODELINE_POINT3]) != 0.0
        error("Nodes with non-zero third coordinates are not supported in this 2D mesh reader.")
      end
      pts[uint64(toks[TOKN_NODELINE_ELNUM])] = pt
      l = readline(is)
    end
    if !beginswith(l, "\$EndNodes") error("End of nodes section not found in mesh file.") end
    pts
  end
end

function is_polytope_or_endmarker(line::ASCIIString)
  if beginswith(line, "\$EndElements")
    true
  else
    const toks = split(line, ' ')
    !is_lower_order_el_type(eltypenum(toks[TOKN_ELLINE_ELTYPE]))
  end
end

function read_through_line_starting(io::IO, line_start::ASCIIString)
  l = readline(io)
  while l != "" && !beginswith(l, line_start)
    l = readline(io)
  end
  l != ""
end

# Read the stream until the line condition function returns true or the stream is exhausted, returning
# either the line passing the condition, or "" if the stream was exhausted, and (in either case)
# the number of lines read including the successful line if any.
function read_until(io::IO, line_cond_fn::Function)
  l = readline(io)
  lines_read = 1
  while l != "" && !line_cond_fn(l)
    l = readline(io)
    lines_read += 1
  end
  l, lines_read
end

# Mesh File Reading Utilities
###############################################



# Gmsh msh format reading support
# ----------------------------------


