
#file to create the run specific parameters files

#getting default parameters values
param_gen_values=loadtxt(param_folder+'param_gen.csv',dtype=object,delimiter=';',usecols=(1,))
param_model_values=loadtxt(param_folder+'param_model.csv',dtype=object,delimiter=';',usecols=(1,))

#reading the run specific parameters list and updating parameters values
with open(preprod_folder+study+'.csv', 'r') as f:
    reader = csv.reader(f,  delimiter=';')
    for row in reader:
        if row[0]==str(run_no):
            for p in list(range(len(table_param_names))):
                if table_param_names[p] in param_gen_names:
                    param_gen_values[param_gen_names==table_param_names[p]]=row[p+2]
                elif table_param_names[p] in param_model_names:
                    param_model_values[param_model_names==table_param_names[p]]=row[p+2]

#saving parameters files in run folder
savetxt(run_folder+'param_gen.csv', column_stack((param_gen_names, param_gen_values)), delimiter=";", fmt="%s")
savetxt(run_folder+'param_model.csv', column_stack((param_model_names, param_model_values)), delimiter=";", fmt="%s")
