# -*- coding: utf-8 -*-
import os 
class CPSOS():
    def mergeFiles(self,inputdir,output_file=None):
        if output_file is None:
            output_file_path = inputdir+".txt"
        file_names = os.listdir(inputdir)
        with open(output_file_path,'wb') as wfh:
            for file_name in file_names:
                print(file_name)
                file_path = os.path.join(inputdir,file_name)
                if file_name == ".DS_Store" or os.path.isdir(file_path):
                    continue
                with open(file_path) as rfh:
                    for one_line in rfh:
                        wfh.write(one_line)
    def googleDriveFilePath(self,google_drive,file_path,ch="\\"):
        file_path=os.path.join(google_drive,file_path)
        return(file_path.replace("\\",ch))