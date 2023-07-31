if to_run!=[0] and to_run!=[] and study!='run_tests4' and study!='run_tests9' and study!='run_tests11': 
    #checking if for all runs with dim>1 there is one, and only one, previous run defined
    prev_run_id=zeros(shape(run_id))
    for run_no in run_id:
        dim=int(table_param_values[run_id==run_no,table_param_names=='dim'])
        stochastic=int(table_param_values[run_id==run_no,table_param_names=='stochastic'])
        kappa1=float(table_param_values[run_id==run_no,table_param_names=='kappa1'])
        kappa2=float(table_param_values[run_id==run_no,table_param_names=='kappa2'])
        pi1=float(table_param_values[run_id==run_no,table_param_names=='pi1'])
        pi2=float(table_param_values[run_id==run_no,table_param_names=='pi2'])
        if dim!=1:
            run_param_values=table_param_values[run_id==run_no,:]
            prev_run_param_values=run_param_values
            is_prev=ones(shape(run_id))
            if stochastic==0:
                prev_run_param_values[:,table_param_names=='dim']=dim-1
            else:
                prev_run_param_values[:,table_param_names=='stochastic']=0          
            for i in list(range(shape(is_prev)[0])):
                if dim==2:
                    if pi1==0 and pi2==0:
                        for ptest in list_prev1:
                            is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names==ptest]==prev_run_param_values[:,table_param_names==ptest])
                    elif pi2==0:
                        prev_run_param_values[:,table_param_names=='dim']=dim
                        for ptest in list_prev1+['pi2']:
                            is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names==ptest]==prev_run_param_values[:,table_param_names==ptest])
                        sort_pi1=sorted(table_param_values[is_prev==1,table_param_names=='pi1'], key=float)
                        pi1_prev=0
                        for j in list(range(shape(sort_pi1)[0])):
                            if sort_pi1[j]<pi1:
                                pi1_prev=sort_pi1[j]
                        is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names=='pi1']==pi1_prev)
                    else:
                        prev_run_param_values[:,table_param_names=='dim']=dim
                        for ptest in list_prev1+['pi1']:
                            is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names==ptest]==prev_run_param_values[:,table_param_names==ptest])
                        sort_pi2=sorted(table_param_values[is_prev==1,table_param_names=='pi2'], key=float)
                        pi2_prev=0
                        for j in list(range(shape(sort_pi2)[0])):
                            if sort_pi2[j]<pi2:
                                pi2_prev=sort_pi2[j]
                        is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names=='pi2']==pi2_prev)
                if dim==3:
                    if kappa1==0 and kappa2==0:
                        for ptest in list_prev2:
                            is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names==ptest]==prev_run_param_values[:,table_param_names==ptest])
                    elif kappa2==0:
                        prev_run_param_values[:,table_param_names=='dim']=dim
                        for ptest in list_prev2+['kappa2']:
                            is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names==ptest]==prev_run_param_values[:,table_param_names==ptest])
                        sort_kappa1=sorted(table_param_values[is_prev==1,table_param_names=='kappa1'], key=float)
                        kappa1_prev=0
                        for j in list(range(shape(sort_kappa1)[0])):
                            if sort_kappa1[j]<kappa1:
                                kappa1_prev=sort_kappa1[j]
                        is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names=='kappa1']==kappa1_prev)
                    else:
                        prev_run_param_values[:,table_param_names=='dim']=dim
                        for ptest in list_prev2+['kappa1']:
                            is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names==ptest]==prev_run_param_values[:,table_param_names==ptest])
                        sort_kappa2=sorted(table_param_values[is_prev==1,table_param_names=='kappa2'], key=float)
                        kappa2_prev=0
                        for j in list(range(shape(sort_kappa2)[0])):
                            if sort_kappa2[j]<kappa2:
                                kappa2_prev=sort_kappa2[j]
                        is_prev[i]=is_prev[i]*(table_param_values[i,table_param_names=='kappa2']==kappa2_prev)            
                #is_prev[i]=(table_param_values[i,:]==prev_run_param_values).all()
            if sum(is_prev)==0:
                print("No previous run defined for run number "+str(run_no))
                well_defined=0
            elif sum(is_prev)>1:
                print("Several previous runs correpond to run number "+str(run_no)+": runs "+str(run_id[is_prev==True]))
                well_defined=0
            else:
                prev_run_id[run_id==run_no]=int(run_id[is_prev==True]) 
                
if to_run!=[0] and to_run!=[] and study=='run_tests4': 
    prev_run_id=zeros(shape(run_id))
    prev_run_id[1:]=run_id[0:-1]
    prev_run_id[11]=run_id[2]
    prev_run_id[19]=run_id[2]
    prev_run_id[27]=run_id[2]
    prev_run_id[35]=run_id[2]
    prev_run_id[43]=run_id[2]

if to_run!=[0] and to_run!=[] and study=='run_tests11': 
    prev_run_id=zeros(shape(run_id))
    prev_run_id[1:]=run_id[0:-1]
    prev_run_id[10]=run_id[2]
    prev_run_id[18]=run_id[2]
    prev_run_id[24]=run_id[2]
    prev_run_id[28]=run_id[2]
    prev_run_id[41]=run_id[33]
    prev_run_id[49]=run_id[33]
    prev_run_id[55]=run_id[33]
    prev_run_id[59]=run_id[33]
