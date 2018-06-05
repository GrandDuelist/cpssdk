from cpspandas import CPSPandas

cpspandas  = CPSPandas()
def testReadDFFromZipFile():
    tele = "a"; teleName = {"a": "telecom","b":"unicom"}[tele]; date={"a":"2013-08-21.zip","b":"2011-07-27.zip"}[tele]
    filePath = r"I:\A-Research\E-Cellphone\raw\{teleName}\{date}".format(teleName=teleName,date=date)
    df = cpspandas.zipToPandas(filePath,startsWith="part")
    print(df.head)



if __name__=='__main__':
    testReadDFFromZipFile()
