import pandas as pd
from zipfile import ZipFile
import os 

class CPSPandas():
    def __init__(self):
        pass
    
    def zipToPandas(self,filePath,startsWith="", endsWith="",contain="",header=None,error_bad_lines=True):
        zip_file = ZipFile(filePath)
        alldf = []
        for textFile in zip_file.infolist():
            fileName = textFile.filename
            baseName = os.path.basename(fileName)
            if baseName.startswith(startsWith) and baseName.endswith(endsWith) and (contain in baseName):
                try:
                    oneFileDF = pd.read_csv(zip_file.open(fileName),header=header,error_bad_lines=error_bad_lines)
                    alldf.append(oneFileDF)
                except Exception as e:
                    print(e)
        df = pd.concat(alldf)
        return(df)
        

            
            
