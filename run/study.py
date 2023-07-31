# -*- coding: utf-8 -*-

#reading the study
#this is where the users has to define the study to use, the runs to compute and enter a few information
with open(preprod_folder+'study.csv', 'r') as f:
    reader = csv.reader(f,  delimiter=';')
    for row in reader:
        if row[0]=='to_run':
            nb_range_runs=len(row)-1
            for l in list(range(1,len(row[1:])+1)):
                if nb_range_runs==len(row)-1 and row[l]=='':
                    nb_range_runs=l-1
            if nb_range_runs%2==0:
                to_run=[]
                for l in list(range(1,nb_range_runs//2+1)):
                    to_run=to_run+list(range(int(row[2*l-1]),int(row[2*l])+1))
            else:
                print("The runs to compute are not well defined, please give a sequence of first runs and last runs")
        elif row[0]=='to_analyze':
            nb_range_runs=len(row)-1
            for l in list(range(1,len(row[1:])+1)):
                if nb_range_runs==len(row)-1 and row[l]=='':
                    nb_range_runs=l-1
            if nb_range_runs%2==0:
                to_analyze=[]
                for l in list(range(1,nb_range_runs//2+1)):
                    to_analyze=to_analyze+list(range(int(row[2*l-1]),int(row[2*l])+1))
            else:
                print("The runs to analyse are not well defined, please give a sequence of first runs and last runs")
        elif row[0]=='study':
            study=row[1]
        else:
            len_param_list_prev=len(row)
            for l in list(range(1,len(row[1:])+1)):
                if len_param_list_prev==len(row) and row[l]=='':
                    len_param_list_prev=l
            exec(row[0]+'=row[1:len_param_list_prev]')
