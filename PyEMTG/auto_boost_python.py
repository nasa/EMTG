header_files = ["SpacecraftOptions","StageOptions","LaunchVehicleOptions","PowerSystemOptions","PropulsionSystemOptions"]

file_handle = open("../src/Python/PyHardware.cpp",'w')

file_handle.write("")

file_handle.write("#define BOOST_PYTHON_STATIC_LIB\n")
file_handle.write("\n")
file_handle.write('#include "boost/python/module.hpp"\n')
file_handle.write('#include "boost/python/class.hpp"\n')
file_handle.write('#include "boost/python/dict.hpp"\n')
file_handle.write('#include "boost/python/import.hpp"\n')
file_handle.write('#include "boost/python/manage_new_object.hpp"\n')
file_handle.write("\n")
file_handle.write('#include "LaunchVehicleOptions.h"\n')
file_handle.write('#include "LaunchVehicle.h"\n')
file_handle.write("\n")
file_handle.write('#include "PowerSystemOptions.h"\n')
file_handle.write('#include "PowerSystem.h"\n')
file_handle.write("\n")
file_handle.write('#include "PropulsionSystemOptions.h"\n')
file_handle.write('#include "ElectricPropulsionSystem.h"\n')
file_handle.write('#include "ChemicalPropulsionSystem.h"\n')
file_handle.write("\n")
file_handle.write('#include "StageOptions.h"\n')
file_handle.write('#include "Stage.h"\n')
file_handle.write("\n")
file_handle.write('#include "SpacecraftOptions.h"\n')
file_handle.write('#include "Spacecraft.h"\n')
file_handle.write("\n")
file_handle.write("namespace EMTG\n")
file_handle.write("{\n")
file_handle.write("     namespace HardwareModels\n")
file_handle.write("     {\n")
file_handle.write("	    template<class T>\n")
file_handle.write("	    inline PyObject * managingPyObject(T *p)\n")
file_handle.write("	    {\n")
file_handle.write("	        return typename boost::python::manage_new_object::apply<T *>::type()(p);\n")
file_handle.write("	    }\n")
file_handle.write("\n")
file_handle.write("	    template<class Copyable>\n")
file_handle.write("	    boost::python::object\n")
file_handle.write("	    generic__deepcopy__(boost::python::object copyable, boost::python::dict memo)\n")
file_handle.write("	    {\n")
file_handle.write('	        boost::python::object copyMod = boost::python::import("copy");\n')
file_handle.write('	        boost::python::object deepcopy = copyMod.attr("deepcopy");\n')
file_handle.write("\n")
file_handle.write("	        Copyable *newCopyable(new Copyable(boost::python::extract<const Copyable&>(copyable)));\n")
file_handle.write("	        boost::python::object result(boost::python::detail::new_reference(managingPyObject(newCopyable)));\n")
file_handle.write("\n")
file_handle.write("	        // HACK: copyableId shall be the same as the result of id(copyable) in Python -\n")
file_handle.write("	        // please tell me that there is a better way! (and which ;-p)\n")
file_handle.write("	        // int copyableId = (int)(copyable.ptr());\n")
file_handle.write("	        // memo[copyableId] = result;\n")
file_handle.write("\n")
file_handle.write('	        boost::python::extract<boost::python::dict>(result.attr("__dict__"))().update(deepcopy(boost::python::extract<boost::python::dict>(copyable.attr("__dict__"))(),memo));\n')
file_handle.write("\n")
file_handle.write("	        return result;\n")
file_handle.write("	    }\n")
file_handle.write("\n")
file_handle.write("        BOOST_PYTHON_MODULE(PyHardware)\n")
file_handle.write("        {\n")

for header in header_files:
        
    header_handle = open("../src/HardwareModels/" + header + ".h")
    
    header_lines = header_handle.readlines()
    
    in_class = False
    
    for line in header_lines:
        
        if in_class:
        
            if line.lstrip(" ").rstrip(" \r\n") == "public:":
                in_public = True
            elif line.lstrip(" ").rstrip(" \r\n") == "private:":
                in_public = False
            elif line.lstrip(" ").rstrip(" \r\n") == "protected:":
                in_public = False
            elif line.rstrip(" \r\n").lstrip(" ") == "};":
                in_class = False
                in_public = False
                

                file_handle.write('                .def("__deepcopy__",&generic__deepcopy__<' + class_name + '>)\n')
                
                for method in methods:
            
                    file_handle.write('                .def("' + method + '",&' + class_name + '::' + method + ')\n')
                    
                for overload_tuple in overload_list:
                    file_handle.write('                .def("' + overload_tuple[0] + '",static_cast<void (' + class_name + '::*)(' + overload_tuple[1] + ')>(&' + class_name + '::' + overload_tuple[0] + '))\n')
                    
                    
                file_handle.write("				;\n\n")
        
            if in_public:
                                
                if not inline:  
                          
                    if class_name not in line and "(" in line and ")" in line and "//" not in line and (";" in line or "inline" in line):
                        
                        if "inline" in line:
                            if "}" in line and "{" in line:
                                inline = False
                            else:
                                inline = True
                            method_name = line.lstrip(" ").split(" ")[2].split("(")[0]

                        else:
                            method_name = line.lstrip(" ").split(" ")[1].split("(")[0]
                        
                        inputs = line.split("(")[1].split(")")[0]
                                                
                        if method_name not in methods and method_name not in overloads:
                    
                            methods.append(method_name)
                            
                            input_list.update({method_name : inputs})
                                                    
                        else:
                            
                            if method_name in methods:
                                methods.remove(method_name)
                            
                            if method_name not in overloads:
                                overloads.append(method_name)
                                
                                new_tuple = (method_name,input_list[method_name])
                            
                                if new_tuple not in overload_list:
                                    overload_list.append(new_tuple)                                
                            
                            new_tuple = (method_name,inputs)
                            
                            if new_tuple not in overload_list:
                                overload_list.append(new_tuple)
                
                elif line.lstrip(" \t").rstrip(" \r\n\t") == "}":
            
                    inline = False
                    
        
                # if "void set" in line:
#
#                     methodName = line.split("set")[1].split("(")[0]
#
#                     file_handle.write('                .def("set' + methodName + '",&' + class_name + '::set' + methodName + ')\n')
#
#                 elif "get" in line and ";" in line:
#
#                     methodName = line.split("get")[1].split("(")[0]
#
#                     file_handle.write('                .def("get' + methodName + '",&' + class_name + '::get' + methodName + ')\n')

        else:
            if "class" in line and "//" not in line:
                class_name = line.split("class")[1].lstrip(" ").rstrip(" \r\n")
                
                file_handle.write('			boost::python::class_<' + class_name + '>("' + class_name + '",boost::python::init<const std::string&>())\n')
                
                in_class = True
                
                in_public = False
                
                inline = False
                
                methods = []
                                
                overloads = []
                
                overload_list = []
                
                input_list = {}
    
file_handle.write("        }\n")
file_handle.write("    }\n")
file_handle.write("}\n")

file_handle.close()