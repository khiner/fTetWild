// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/CSGTreeParser.hpp>

#include <floattetwild/Timer.h>
#include <floattetwild/Logger.hpp>

#include <floattetwild/geo_mesh_reorder.h>

#include <cstring>
#include <iterator>

namespace floatTetWild {

namespace {
// Build a geo::Mesh from vertex and face vectors, carry the tags through as a facet
// attribute, Morton reorder it, and read the permuted result back out. Named load_mesh
// on MeshIO, but it touches no files.
void build_reordered_sf_mesh(std::vector<Vector3>&  points,
                             std::vector<Vector3i>& faces,
                             geo::Mesh&             input,
                             std::vector<int>&      flags)
{
    logger().debug("Loading mesh from data...");
    Timer timer;
    timer.start();
    input.clear(false, false);

    input.clear();
    input.vertices.create_vertices((int)points.size());
    for (int i = 0; i < (int)input.vertices.nb(); ++i) {
        geo::vec3& p = input.vertices.point(i);
        p[0]         = points[i](0);
        p[1]         = points[i](1);
        p[2]         = points[i](2);
    }
    // Setup faces
    input.facets.create_triangles((int)faces.size());

    for (int c = 0; c < (int)input.facets.nb(); ++c) {
        for (int lv = 0; lv < 3; ++lv) {
            input.facets.set_vertex(c, lv, faces[c](lv));
        }
    }

    bool is_valid = (flags.size() == input.facets.nb());
    if (is_valid) {
        assert(flags.size() == input.facets.nb());
        geo::Attribute<int> bflags(input.facets.attributes(), "bbflags");
        for (int index = 0; index < (int)input.facets.nb(); ++index) {
            bflags[index] = flags[index];
        }
    }

    geo::mesh_reorder(input, geo::MESH_ORDER_MORTON);

    if (is_valid) {
        flags.clear();
        flags.resize(input.facets.nb());
        geo::Attribute<int> bflags(input.facets.attributes(), "bbflags");
        for (int index = 0; index < (int)input.facets.nb(); ++index) {
            flags[index] = bflags[index];
        }
    }

    points.resize(input.vertices.nb());
    for (size_t i = 0; i < points.size(); i++)
        points[i] << (input.vertices.point(i))[0], (input.vertices.point(i))[1],
          (input.vertices.point(i))[2];

    faces.resize(input.facets.nb());
    for (size_t i = 0; i < faces.size(); i++)
        faces[i] << int(input.facets.vertex(i, 0)), int(input.facets.vertex(i, 1)),
          int(input.facets.vertex(i, 2));
}
// Just enough JSON for the csg tree format. Objects and strings are read straight into CSGTree, and
// a value under any other key is skipped whatever its type, which is what parsing the file into a
// general json value and reading three keys back out of it used to do.
class CsgReader
{
  public:
    explicit CsgReader(std::istream& in)
        : text_((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>())
    {}

    bool parse(CSGTree& tree)
    {
        if (!parse_node(tree))
            return false;
        skip_space();
        if (pos_ != text_.size())
            return fail("trailing content after the tree");
        return true;
    }

    const std::string& error() const { return error_; }

  private:
    bool fail(const std::string& what)
    {
        if (error_.empty())
            error_ = what + ", at character " + std::to_string(pos_ + 1);
        return false;
    }

    void skip_space()
    {
        while (pos_ < text_.size() && (text_[pos_] == ' ' || text_[pos_] == '\t' ||
                                       text_[pos_] == '\n' || text_[pos_] == '\r'))
            ++pos_;
    }

    bool at(char c)
    {
        skip_space();
        return pos_ < text_.size() && text_[pos_] == c;
    }

    bool take(char c)
    {
        if (!at(c))
            return false;
        ++pos_;
        return true;
    }

    // The four hex digits of a \u escape, joined with a following low surrogate when there is one,
    // written out as UTF-8.
    bool parse_unicode_escape(std::string& out)
    {
        unsigned int code = 0;
        if (!parse_hex4(code))
            return false;
        if (code >= 0xD800 && code <= 0xDBFF) {
            if (pos_ + 1 < text_.size() && text_[pos_] == '\\' && text_[pos_ + 1] == 'u') {
                pos_ += 2;
                unsigned int low = 0;
                if (!parse_hex4(low))
                    return false;
                if (low < 0xDC00 || low > 0xDFFF)
                    return fail("a high surrogate is not followed by a low one");
                code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
            }
            else
                return fail("a high surrogate is not followed by a low one");
        }

        if (code < 0x80)
            out.push_back(char(code));
        else if (code < 0x800) {
            out.push_back(char(0xC0 | (code >> 6)));
            out.push_back(char(0x80 | (code & 0x3F)));
        }
        else if (code < 0x10000) {
            out.push_back(char(0xE0 | (code >> 12)));
            out.push_back(char(0x80 | ((code >> 6) & 0x3F)));
            out.push_back(char(0x80 | (code & 0x3F)));
        }
        else {
            out.push_back(char(0xF0 | (code >> 18)));
            out.push_back(char(0x80 | ((code >> 12) & 0x3F)));
            out.push_back(char(0x80 | ((code >> 6) & 0x3F)));
            out.push_back(char(0x80 | (code & 0x3F)));
        }
        return true;
    }

    bool parse_hex4(unsigned int& code)
    {
        if (pos_ + 4 > text_.size())
            return fail("a \\u escape is cut short");
        code = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = text_[pos_++];
            code <<= 4;
            if (c >= '0' && c <= '9')
                code |= unsigned(c - '0');
            else if (c >= 'a' && c <= 'f')
                code |= unsigned(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F')
                code |= unsigned(c - 'A' + 10);
            else
                return fail("a \\u escape has a non-hex digit");
        }
        return true;
    }

    bool parse_string(std::string& out)
    {
        if (!take('"'))
            return fail("expected a string");
        out.clear();
        while (pos_ < text_.size()) {
            const char c = text_[pos_++];
            if (c == '"')
                return true;
            if (c != '\\') {
                out.push_back(c);
                continue;
            }
            if (pos_ >= text_.size())
                break;
            switch (text_[pos_++]) {
            case '"':
                out.push_back('"');
                break;
            case '\\':
                out.push_back('\\');
                break;
            case '/':
                out.push_back('/');
                break;
            case 'b':
                out.push_back('\b');
                break;
            case 'f':
                out.push_back('\f');
                break;
            case 'n':
                out.push_back('\n');
                break;
            case 'r':
                out.push_back('\r');
                break;
            case 't':
                out.push_back('\t');
                break;
            case 'u':
                if (!parse_unicode_escape(out))
                    return false;
                break;
            default:
                return fail("unknown escape in a string");
            }
        }
        return fail("unterminated string");
    }

    // Consumes any value, for keys the tree does not use.
    bool skip_value()
    {
        skip_space();
        if (pos_ >= text_.size())
            return fail("expected a value");
        if (at('"')) {
            std::string ignored;
            return parse_string(ignored);
        }
        if (at('{') || at('[')) {
            const bool is_object = at('{');
            const char close     = is_object ? '}' : ']';
            ++pos_;
            if (take(close))
                return true;
            do {
                if (is_object) {
                    std::string key;
                    if (!parse_string(key))
                        return false;
                    if (!take(':'))
                        return fail("expected ':' after a key");
                }
                if (!skip_value())
                    return false;
            } while (take(','));
            return take(close) ? true : fail("unterminated object or array");
        }
        // A number, true, false or null. None of them appear in the tree itself, so the exact
        // spelling does not matter, only where it ends.
        const size_t start = pos_;
        while (pos_ < text_.size() && std::strchr(",}] \t\n\r", text_[pos_]) == nullptr)
            ++pos_;
        return pos_ > start ? true : fail("expected a value");
    }

    bool parse_child(CSGTree::Child& child)
    {
        if (at('"')) {
            child.kind = CSGTree::Child::Mesh;
            return parse_string(child.mesh);
        }
        child.kind    = CSGTree::Child::Subtree;
        child.subtree = std::make_shared<CSGTree>();
        return parse_node(*child.subtree);
    }

    bool parse_node(CSGTree& tree)
    {
        if (!take('{'))
            return fail("expected an object");

        bool have_operation = false, have_left = false, have_right = false;
        if (!take('}')) {
            do {
                std::string key;
                if (!parse_string(key))
                    return false;
                if (!take(':'))
                    return fail("expected ':' after a key");

                if (key == "operation") {
                    if (!parse_string(tree.operation))
                        return false;
                    have_operation = true;
                }
                else if (key == "left") {
                    if (!parse_child(tree.left))
                        return false;
                    have_left = true;
                }
                else if (key == "right") {
                    if (!parse_child(tree.right))
                        return false;
                    have_right = true;
                }
                else if (!skip_value())
                    return false;
            } while (take(','));

            if (!take('}'))
                return fail("expected ',' or '}'");
        }

        if (!have_operation || !have_left || !have_right)
            return fail("a node needs \"operation\", \"left\" and \"right\"");
        if (tree.operation != "union" && tree.operation != "intersection" &&
            tree.operation != "difference")
            return fail("unknown operation \"" + tree.operation + "\"");
        return true;
    }

    const std::string text_;
    size_t            pos_ = 0;
    std::string       error_;
};

}  // namespace

bool CSGTreeParser::parse(std::istream& in, CSGTree& tree, std::string& error)
{
    CsgReader reader(in);
    if (reader.parse(tree))
        return true;
    error = reader.error();
    return false;
}

void CSGTreeParser::get_meshes_aux(const CSGTree&              csg_tree_node,
                                   std::vector<std::string>&   meshes,
                                   std::map<std::string, int>& existings,
                                   int&                        index,
                                   CSGTree&                    current_node)
{
    // A name already seen keeps the index it was given the first time, so a mesh used twice is
    // only loaded once.
    const auto resolve = [&](const CSGTree::Child& child) {
        CSGTree::Child out;
        if (child.kind == CSGTree::Child::Subtree) {
            out.kind    = CSGTree::Child::Subtree;
            out.subtree = std::make_shared<CSGTree>();
            get_meshes_aux(*child.subtree, meshes, existings, index, *out.subtree);
            return out;
        }

        const auto iter = existings.find(child.mesh);
        if (iter == existings.end()) {
            meshes.push_back(child.mesh);
            existings[child.mesh] = index;
            out.id                = index;
            ++index;
        }
        else
            out.id = iter->second;

        out.kind = CSGTree::Child::Id;
        return out;
    };

    current_node.operation = csg_tree_node.operation;
    current_node.left      = resolve(csg_tree_node.left);
    current_node.right     = resolve(csg_tree_node.right);
}

void CSGTreeParser::merge(const std::vector<std::vector<Vector3>>&  Vs,
                          const std::vector<std::vector<Vector3i>>& Fs,
                          std::vector<Vector3>&                     V,
                          std::vector<Vector3i>&                    F,
                          geo::Mesh&                                sf_mesh,
                          std::vector<int>&                         tags)
{
    V.clear();
    F.clear();

    for (const auto& vv : Vs)
        V.insert(V.end(), vv.begin(), vv.end());

    int offset = 0;
    int size   = 0;
    for (int id = 0; id < Fs.size(); ++id) {
        const auto& ff = Fs[id];

        for (const auto fid : ff)
            F.push_back(Vector3i(fid(0) + offset, fid(1) + offset, fid(2) + offset));

        tags.insert(tags.begin() + size, Fs[id].size(), id);
        size += Fs[id].size();
        offset += Vs[id].size();
    }

    build_reordered_sf_mesh(V, F, sf_mesh, tags);
}

void CSGTreeParser::get_max_id_aux(const CSGTree& csg_tree_node, int& max)
{
    for (const CSGTree::Child* child : {&csg_tree_node.left, &csg_tree_node.right}) {
        if (child->kind == CSGTree::Child::Id)
            max = std::max(max, child->id);
        else
            get_max_id_aux(*child->subtree, max);
    }
}

bool CSGTreeParser::keep_tet(const CSGTree&                      csg_tree_with_ids,
                             const int                           t_id,
                             const std::vector<VectorXd>&        w)
{
    const auto inside = [&](const CSGTree::Child& child) {
        return child.kind == CSGTree::Child::Id ? w[child.id][t_id] > 0.5
                                                : keep_tet(*child.subtree, t_id, w);
    };

    const bool left_inside  = inside(csg_tree_with_ids.left);
    const bool right_inside = inside(csg_tree_with_ids.right);

    const std::string& op = csg_tree_with_ids.operation;
    if (op == "union")
        return left_inside || right_inside;
    if (op == "intersection")
        return left_inside && right_inside;
    if (op == "difference")
        return left_inside && !right_inside;

    assert(false);
    return false;
}

}  // namespace floatTetWild
