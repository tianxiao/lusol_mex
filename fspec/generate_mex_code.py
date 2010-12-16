#!/usr/bin/env python

def get_c_symb(symb):
    a = symb + "_FUNC"
    return a.upper()

def get_f_symb(symb):
    return symb + "_"

def write_subroutine_call(c_file,fspec):
    spec_filename = fspec + ".fspec"
    spec_file = open(spec_filename,"r")

    base_error = "lusol_mex"

    spec_lines = list()

    for line in spec_file:
        spec_lines.append(line.split())

    func_name = spec_lines[0][0]

    spec_lines.pop(0)
    spec_lines.pop(0)

    num_inputs = len(spec_lines)

    # start writing the file

    c_file.write("void gateway_" + func_name + "(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {\n")

    # check number of right hand side arguments
    c_file.write("  /* check nrhs */\n")
    c_file.write("  if (nrhs != {0}) {{\n".format(num_inputs+1))
    c_file.write("    mexErrMsgIdAndTxt(\"{0}:{1}\",\"incorrect number of inputs for {1}.\");\n".format(base_error,func_name))
    c_file.write("  }\n\n")

    var_counter = 1
    var_list = list()
    for var_data in spec_lines:
        var_name = var_data[0]
        var_list.append(var_name)
        c_type = var_data[1]
        m_type = var_data[2]

        error_stmt = "{0} must have type {1} and size at least ({2},{3})".format(var_name,m_type,var_data[3],var_data[4])

        # modify M and N
        try:
            M = int(var_data[3])
        except ValueError:
            M = "*" + var_data[3]
        try:
            N = int(var_data[4])
        except ValueError:
            N = "*" + var_data[4]

        c_file.write("  /* check prhs[{0}], variable: {1} */\n".format(var_counter,var_name))
        # construct code
        out_str = ""
        out_str = out_str + "  if ("
        out_str = out_str + "!mxIsClass(prhs[{0}],\"{1}\")".format(var_counter,m_type)
        out_str = out_str + " || "
        out_str = out_str + "mxGetM(prhs[{0}]) < {1}".format(var_counter,M)
        out_str = out_str + " || "
        out_str = out_str + "mxGetN(prhs[{0}]) != {1}".format(var_counter,N)
        out_str = out_str + ") { \n" # end of if statement
        out_str = out_str + "    mexErrMsgIdAndTxt(\"{0}:{1}\",\"{2}\");\n".format(base_error,func_name,error_stmt)
        out_str = out_str + "  }\n" # close the if block

        out_str = out_str + "  {0} *{1} = ({0}*) mxGetData(prhs[{2}]);\n\n".format(c_type,var_name,var_counter)
        c_file.write(out_str)
        var_counter = var_counter + 1


    # write the function call

    arg_list = ""
    for var_data in spec_lines:
        arg_list = arg_list + var_data[0] + ", "
    arg_list = arg_list[:-2]
    c_file.write("  {0}({1});\n".format(get_c_symb(fspec),arg_list))
    c_file.write("}\n\n")

    # close the spec file
    spec_file.close()


def declare_subroutine_call(c_file,fspec):
    spec_filename = fspec + ".fspec"
    spec_file = open(spec_filename,"r")

    spec_lines = list()
    for line in spec_file:
        spec_lines.append(line.split())

    func_name = spec_lines[0][0]
    spec_lines.pop(0)
    spec_lines.pop(0)

    arg_str = ""

    for var_data in spec_lines:
        arg_str = arg_str + "{0} *{1}, ".format(var_data[1],var_data[0])
    arg_str = arg_str[:-2]
    
    c_file.write("void {0}({1});\n".format(get_c_symb(func_name),arg_str))

    # close the spec file
    spec_file.close()


# output the c code

fspec_list_filename = "fspec_list.txt"
fspec_list_file = open(fspec_list_filename,"r")
fspec_list = fspec_list_file.readlines()

c_filename = "../lusol_mex.c"
c_file = open(c_filename,"w")

# start the c file
c_file.write("#include \"mex.h\"\n")
c_file.write("#include <string.h>\n")
c_file.write("#include <stdlib.h>")
c_file.write("\n\n")

# generate define statements. I am doing this in case some other fortran compiler
# generates different symbols for fortran subroutines
for func in fspec_list:
    func = func[:-1]
    c_symb = get_c_symb(func)
    f_symb = get_f_symb(func)
    c_file.write("#define {0} {1}\n".format(c_symb,f_symb))

# declare the functions
c_file.write("\n")
for func in fspec_list:
    gateway_func = "gateway_" + func[:-1]
    c_file.write("void {0}(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);\n".format(gateway_func))

# declare the subroutine calls
c_file.write("\n")
for func in fspec_list:
    func = func[:-1]
    declare_subroutine_call(c_file,func)

# output the main gateway function
c_file.write("\n")
c_file.write("void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { \n")

# check to make sure the first argument is a string
c_file.write("  if (nrhs < 1) {\n")
c_file.write("    mexErrMsgIdAndTxt(\"lusol_mex:gateway\",\"there are no input arguments.\");\n")
c_file.write("  }\n\n")
c_file.write("  if (!mxIsChar(prhs[0])) {\n")
c_file.write("    mexErrMsgIdAndTxt(\"lusol_mex:gateway\",\"the first argument must be a string indicating the subroutine to call.\");\n")
c_file.write("  }\n\n")

c_file.write("  mwSize strM = mxGetM(prhs[0]);\n")
c_file.write("  mwSize strN = mxGetN(prhs[0]);\n")

c_file.write("  if ( strM != 1 || strN < 1 ) {\n")
c_file.write("    mexErrMsgIdAndTxt(\"lusol_mex:gateway\",\"the input string has incorrect size.\");\n")
c_file.write("  }\n")

c_file.write("  char *subroutine_name = (char*) mxMalloc(sizeof(char)*(strN+1));\n")
c_file.write("  if (!subroutine_name) mexErrMsgIdAndTxt(\"lusol_mex:gateway\",\"could not allocate memory for string.\");\n")
c_file.write("  int check;\n")
c_file.write("  check = mxGetString(prhs[0],subroutine_name,strN+1);\n")
c_file.write("  if (check) {\n")
c_file.write("    mexErrMsgIdAndTxt(\"lusol_mex:gateway\",\"could not read input string.\");\n")
c_file.write("  }\n\n")

# output calls to sub functions
for func in fspec_list:
    func = func[:-1]
    gateway_func = "gateway_" + func
    c_file.write("  if (strcmp(subroutine_name,\"{0}\") == 0) {{\n".format(func))
    c_file.write("    {0}(nlhs, plhs, nrhs, prhs);\n".format(gateway_func))
    c_file.write("    mxFree(subroutine_name);\n")
    c_file.write("    return;\n")
    c_file.write("  }\n")


c_file.write("  mxFree(subroutine_name);\n")
c_file.write("  mexErrMsgIdAndTxt(\"lusol_mex:gateway\",\"the input string does not match any subroutine.\");\n")

c_file.write("}\n")

# output direct calls to subroutines
c_file.write("\n")
for func in fspec_list:
    func = func[:-1]
    write_subroutine_call(c_file,func)
    
