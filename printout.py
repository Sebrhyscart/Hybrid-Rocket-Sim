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

output_file_name = "verbose.o"
data_file_name = "data.o"

def print_title():
    global header
    print(header)

def new_file():
    global output_file_name
    """
    Removes the old instance of the file with the given name, if it exists, 
    so that a new file can be created.
    """
    try:
        if os.path.exists(output_file_name):
            os.remove(output_file_name)
            print(f"Existing file '{output_file_name}' removed.")
        else:
            print(f"File '{output_file_name}' does not exist, ready to create a new one.")
    except Exception as e:
        print(f"An error occurred while trying to remove the file: {e}")

def new_data():
    global data_file_name
    """
    Removes the old instance of the file with the given name, if it exists, 
    so that a new file can be created.
    """
    try:
        if os.path.exists(data_file_name):
            os.remove(data_file_name)
            print(f"Existing file '{data_file_name}' removed.")
        else:
            print(f"File '{data_file_name}' does not exist, ready to create a new one.")
    except Exception as e:
        print(f"An error occurred while trying to remove the file: {e}")

def print_file(*args, sep=' ', end='\n'):
    """
    Mimics the behavior of the print function but writes to a global file.

    Parameters:
    *args: Variable arguments to be written to the file.
    sep (str): Separator to be used between arguments. Default is a single space.
    end (str): String appended after the last argument. Default is a newline. 
    """
    global output_file_name
    try:
        # Open the file in append mode
        with open(output_file_name, 'a') as file:
            # Join the arguments with the separator and append the end string
            file.write(sep.join(map(str, args)) + end)
    except Exception as e:
        print(f"An error occurred in writing output file: {e}")

def print_data(*args, sep=' ', end='\n'):
    """
    Mimics the behavior of the print function but writes to a global file.

    Parameters:
    *args: Variable arguments to be written to the file.
    sep (str): Separator to be used between arguments. Default is a single space.
    end (str): String appended after the last argument. Default is a newline. 
    """
    global data_file_name
    try:
        # Open the file in append mode
        with open(data_file_name, 'a') as file:
            # Join the arguments with the separator and append the end string
            file.write(sep.join(map(str, args)) + end)
    except Exception as e:
        print(f"An error occurred in writing output file: {e}")
