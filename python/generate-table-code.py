import build_cut_table as tt
import argparse

import os


def parse_args():
    parser = argparse.ArgumentParser(description='Export table as a c file')
    parser.add_argument('--filename', type=str, help='name of the file')

    return parser.parse_args()


def emit_table(class_name, ret, fun_name, rows, format_conf, assert_non_empty=True):
    """One `static const std::array<ret, 64> table` and the accessor that indexes it.

    rows[i] is None where the configuration does not exist, and otherwise the list of
    configurations for index i; format_conf turns one of those into its initialiser list.
    """
    cpp = "\tconst {}& {}::{} {{\n".format(ret, class_name, fun_name)
    cpp += "\t\tstatic const std::array<{}, {}> table= {{{{\n".format(ret, len(rows))
    for confs in rows:
        cpp += "\n"
        if confs is None:
            cpp += "\t\t\t{},"
            continue

        cpp += "\t\t\t{\n"
        for c in confs:
            body = format_conf(c)
            # Each entry is written with a trailing comma and the last one is chopped back off,
            # which for an empty configuration eats one of the indent tabs instead.
            cpp += "\n\t\t\t\t{\n\t\t\t\t\t" + (body + "," if body else "")
            cpp = cpp[:-1]
            cpp += "\n\t\t\t\t},"
        cpp = cpp[:-1]
        cpp += "\n\t\t\t},"
    cpp = cpp[:-1]
    cpp += "\n\t\t}};\n\n"

    cpp += "\t\t{}assert(!table[idx].empty());\n".format("" if assert_non_empty else "//")
    cpp += "\t\treturn table[idx];\n"
    cpp += "\t}\n\n\n"
    return cpp


def format_tets(c):
    # A configuration that reaches vertex 9 without using 8 is renumbered down onto 8.
    max_v = max(max(tet) for tet in c)
    has_8 = any(8 in tet for tet in c)
    decrease_9 = max_v == 9 and not has_8
    fix = lambda v: 8 if decrease_9 and v == 9 else v
    return ",".join("Vector4i({}, {}, {}, {})".format(*[fix(v) for v in tet]) for tet in c)


def format_diags(c):
    return ",".join("Vector2i({}, {})".format(e[0], e[1]) for e in c)


def format_surface(c):
    return ",".join(
        "{{{{{}, {}, {}, {}}}}}".format(*[str(b).lower() for b in t]) for t in c)


def format_face_ids(c):
    return ",".join("Vector4i({}, {}, {}, {})".format(*f) for f in c)


def main():
    class_name = "CutTable"

    args = parse_args()
    save_file = None
    if args.filename:
        fname = os.path.basename(args.filename)
        save_file = os.path.join(os.getcwd(), args.filename)
    else:
        fname = "test"

    table = tt.CutTable()

    # (return type, accessor, data, formatter, assert that the entry is non-empty)
    tables = [
        ("std::vector<std::vector<Vector4i>>", "get_tet_confs(const int idx)",
         table.table, format_tets, True),
        ("std::vector<std::vector<Vector2i>>", "get_diag_confs(const int idx)",
         table.edges_table, format_diags, False),
        ("std::vector<std::vector<std::array<bool, 4>>>", "get_surface_conf(const int idx)",
         table.track_faces, format_surface, True),
        ("std::vector<std::vector<Vector4i>>", "get_face_id_conf(const int idx)",
         table.original_faces, format_face_ids, True),
    ]

    hpp = "#pragma once\n\n"
    hpp += "#include <floattetwild/Types.hpp>\n\n"
    hpp += "#include <array>\n"
    hpp += "#include <vector>\n\n"

    hpp += "namespace floatTetWild {{\n\tclass {} {{\npublic:\n".format(class_name)
    for ret, fun_name, _, _, _ in tables:
        hpp += "\t\tstatic const {}& {};\n".format(ret, fun_name)

    # The per-configuration accessors, which just index the ones above.
    for elem, fun_name, body in [
        ("std::vector<Vector4i>", "get_tet_conf(const int idx, const int cfg)", "get_tet_confs"),
        ("std::vector<std::array<bool, 4>>", "get_surface_conf(const int idx, const int cfg)",
         "get_surface_conf"),
        ("std::vector<Vector4i>", "get_face_id_conf(const int idx, const int cfg)",
         "get_face_id_conf"),
    ]:
        hpp += "\t\tstatic inline const {}& {}{{ return {}(idx)[cfg]; }}\n".format(
            elem, fun_name, body)

    hpp += "\t};\n}\n"

    cpp = '#include "{}.hpp"\n\n'.format(fname)
    cpp += '#include <cassert>\n\n'

    cpp += "namespace floatTetWild {\n"
    for ret, fun_name, data, format_conf, assert_non_empty in tables:
        # Which entries exist is decided by the tet-configuration table throughout.
        rows = [None if table.table[i] is None else data[i] for i in range(len(data))]
        cpp += emit_table(class_name, ret, fun_name, rows, format_conf, assert_non_empty)
    cpp += "}\n"

    if save_file:
        with open(save_file + ".hpp", 'w') as f:
            f.write(hpp)

        with open(save_file + ".cpp", 'w') as f:
            f.write(cpp)
    else:
        print(cpp)


if __name__ == '__main__':
    main()
