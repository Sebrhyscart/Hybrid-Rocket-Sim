import os

# ASCII art header
header = '''               .-==:                           
              :=:  :-                 __  __        __  __              _           
             :=:    .-               |  \/  |      |  \/  |            | |    
             -:      .:              | \  / |  ___ | \  / |  __ _  ___ | |_  ___  _ __ 
            :-.       .              | |\/| | / __|| |\/| | / _` |/ __|| __|/ _ \| '__| 
            -:         .             | |  | || (__ | |  | || (_| |\__ \| |_|  __/| |   
           .==-        .             |_|__|_| \___||_|  |_| \__,_||___/ \__|\___||_|   
           :=++=:                    |  __ \             | |        | |                 
           :-==:=-                   | |__) | ___    ___ | | __ ___ | |_  _ __  _   _ 
           -::=: :=:                 |  _  / / _ \  / __|| |/ // _ \| __|| '__|| | | | 
           :: -:   ::.               | | \ \| (_) || (__ |   <|  __/| |_ | |   | |_| | 
           .  :=     ..              |_|  \_\\\\___/  \___||_|\_\\\\___| \__||_|    \__, |
          ..   =.      ..                                                        __/ |
         :+*-  ::     .+*-                                                      |___/ 
        :*#%#=. :    :*%##=.        
       -*%#%##+::.  :###%##*:       =====================================================
     .=##*-:=##*:. =##*-:=##*:           Hybrid Rocket Engine Combustion Simulation
    .+#*:    .=*#=+#*:    .=*#=     
   :**:        .+%#-        .-*+.   
  :+:          -=:-+:          -=:  
 ::           :.    :.           :. 
 ========================================================================================'''

def print_title():
    global header
    print(header)

def new_file(filename):
    """
    Removes the old instance of the file with the given name, if it exists, 
    so that a new file can be created.

    Parameters:
    filename: Name of the file to write over
    """
    if filename != None:
        output_name = "verbose_" + filename + ".o"
    else:
        output_name = "verbose.o"

    try:
        if os.path.exists(output_name):
            os.remove(output_name)
    except Exception as e:
        print(f"An error occurred while trying to remove the file: {e}")

def new_data(filename):
    """
    Removes the old instance of the file with the given name, if it exists, 
    so that a new file can be created.

    Parameters:
    filename: Name of the file to write over
    """
    if filename != None:
        output_name = filename + ".o"
    else:
        output_name = "data.o"

    try:
        if os.path.exists(output_name):
            os.remove(output_name)
    except Exception as e:
        print(f"An error occurred while trying to remove the file: {e}")

def print_file(filename, *args, sep=' ', end='\n'):
    """
    Mimics the behavior of the print function but writes to a global file.

    Parameters:
    filename: Name of the file to write to
    *args: Variable arguments to be written to the file.
    sep (str): Separator to be used between arguments. Default is a single space.
    end (str): String appended after the last argument. Default is a newline. 
    """
    if filename != None:
        output_name = "verbose_printout_" + filename + ".o"
    else:
        output_name = "verbose_printout.o"

    try:
        # Open the file in append mode
        with open(output_name, 'a') as file:
            # Join the arguments with the separator and append the end string
            file.write(sep.join(map(str, args)) + end)
    except Exception as e:
        print(f"An error occurred in writing output file: {e}")

def print_data(filename, *args, sep=' ', end='\n'):
    """
    Mimics the behavior of the print function but writes to a global file.

    Parameters:
    filename: Name of the file to write to
    *args: Variable arguments to be written to the file.
    sep (str): Separator to be used between arguments. Default is a single space.
    end (str): String appended after the last argument. Default is a newline. 
    """
    if filename != None:
        output_name = filename + ".o"
    else:
        output_name = "data.o"

    try:
        # Open the file in append mode
        with open(output_name, 'a') as file:
            # Join the arguments with the separator and append the end string
            file.write(sep.join(map(str, args)) + end)
    except Exception as e:
        print(f"An error occurred in writing output file: {e}")
