#!/home/guorc/bin/python3
# -*- coding: UTF-8 -*-

import argparse,sys,os

def code_help():
        
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help="input table1, required")
    parser.add_argument("-b", help="input table2, required")
    parser.add_argument("-o", help="input name of outfile, required")
    parser.add_argument("-m", default='table_row',choices=['table_column','table_row','table_merge'], help="pick out the needed rows or columns according to table2, default table_row")
    parser.add_argument("-r", action="store_true", help="whether to output variables that aren't in table2, default closed")
    parser.add_argument("-w", default='left', help="merge type [left/right/inner/outer], default left")
    parser.add_argument("-x", type=int, default=0, help="column for table1, default 0")
    parser.add_argument("-y", type=int, default=0, help="column for table2, default 0")
    parser.add_argument("-i", default=None, help="whether or not table1 has header, default None, e.g. 0,1,2 ...")
    parser.add_argument("-j", default=None, help="whether or not table2 has header, default None, e.g. 0,1,2 ...")
    parser.add_argument("-f", default='NA', help="replace this of vacant position, default NA")
    
    args = parser.parse_args()
    script_path = os.path.abspath(sys.argv[0])
    
    if not args.a or not args.b or not args.o:
        print("\nMissing parameters\n")
        os.system('python3 {SP} -h'.format(SP=script_path))
        sys.exit()    

    
    return (args.a,args.b,args.o,args.x,args.y,args.m,args.r,args.i, args.j, args.w,args.f)

aTable1, aTable2,outfile, aTcolumn, aLcolumn, aMethod, aReverse, arHead, alHead, aMergehow, aFill = code_help()

#=========================================================================================================================
import pandas as pd


def opentab(Table,Header,Sep='\t'):
    Table = pd.read_csv(Table,header=Header,sep=Sep)
    return(Table)

    
def Row(Table,Map,Tab_num,Map_num,Rev=False):
    if not Rev:
        Table = Table[Table.iloc[:,Tab_num].isin(Map.iloc[:,Map_num])]
    else:
        Table = Table[~Table.iloc[:,Tab_num].isin(Map.iloc[:,Map_num])]
    return(Table)


def Col(Table,Map,Tab_num,Map_num,Rev=False):
    lsts = [x for x in Map.iloc[:,Map_num] if x in Table.columns]
    
    if Map.shape[0]!=len(lsts):
        print("There are some samples that aren't in table1!")
    
    if not Rev:
        lsts.insert(0,Table.columns[Tab_num])
    else:
        lsts = [x for x in Table.columns if x not in lsts]
    Table = Table.loc[:,lsts]
    return(Table)


def Merge(Table,Map,Tab_num,Map_num,Merge_how='left'):
    Table_names = list(Table.columns)
    Map_names = list(Map.columns)

    Table.set_index(Table_names[Tab_num],inplace=True)
    Map.set_index(Map_names[Map_num],inplace=True)
    del Map_names[Map_num]
    ff = Table_names[Tab_num]
    del Table_names[Tab_num]
    Table = Table.join(Map, how=Merge_how, sort=False, lsuffix='_left', rsuffix='_right')
    aa = []
    for x in Table_names:
        if x in Map_names:
            aa.append(str(x)+'_left')
        else:
            aa.append(x)
    for x in Map_names:
        if x in Table_names:
            aa.append(str(x)+'_right')
        else:
            aa.append(x)
    Table = Table[aa]
    Table.columns = [Table_names+Map_names]
    Table.insert(Tab_num,ff,Table.index,allow_duplicates=True)
    #Table.fillna(Fill,inplace=True)
    return(Table)

def list_merge(Tab,lst,Tab_num,Lst_num,How):
    lst.rename(columns = {x:str(x) for x in lst.columns},inplace=True)
    Tab.rename(columns = {x:str(x) for x in Tab.columns},inplace=True)
    lst.rename(columns = {lst.columns[Lst_num]:'LstID'},inplace=True)
    Torg = Tab.columns
    Tab.rename(columns = {x:"Tab_"+x for x in Torg},inplace=True)
    Tab.rename(columns = {Tab.columns[Tab_num]:'LstID'},inplace=True)
    Tab = pd.merge(Tab,lst,how=How,on='LstID',sort=False)
    Tab.rename(columns = {'LstID':"Tab_"+Torg[Tab_num]},inplace=True)
    TT = Tab.loc[:,["Tab_"+ x for x in Torg]]
    f = [x for x in Tab.columns if x not in TT.columns]
    for x in f:
        TT[x]=Tab.loc[:,x]
    TT.rename(columns = {"Tab_"+x:x for x in Torg},inplace=True)
    #TT.fillna(Fill,inplace=True)
    return(TT)
#=========================================================================================================================

if arHead is None:
    tab1 = opentab(aTable1,arHead)
else:
    tab1 = opentab(aTable1,int(arHead))

if alHead is None:
    tab2 = opentab(aTable2,alHead)
else:
    tab2 = opentab(aTable2,int(alHead))


if aMethod=='table_row':
    tab = Row(tab1,tab2,aTcolumn,aLcolumn,aReverse)
elif aMethod=='table_column': 
    tab = Col(tab1,tab2,aTcolumn,aLcolumn,aReverse)
else:
    tab = Merge(tab1,tab2,aTcolumn,aLcolumn,aMergehow)
    #tab = list_merge(tab1,tab2,aTcolumn,aLcolumn,aMergehow)


if arHead is None:
    tab.to_csv(outfile,sep='\t',index=False,header=0, na_rep=aFill)
else:
    tab.to_csv(outfile,sep='\t',index=False,na_rep=aFill)


