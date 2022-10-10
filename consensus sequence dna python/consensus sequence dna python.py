
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#amin hashemian
import collections
import math
import numpy as np


#dna matrix tanımı
A = np.array([['A', '-', 'T','T', 'G', '-','C', 'T', 'T'], 
              ['A', 'C', 'C','T', 'G', 'G','C', '-', 'T'], 
              ['A', 'C', 'T','A', '-', 'G','C', 'T', 'A'], 
              ['A', 'G', 'T','A', 'G', 'C','A', 'T', 'T']] ) 


#satırlardan elde edilen sonuçları kullanmak içi dizi oluşumu
SOP_Score=[0,0,0,0,0,0,0,0,0,0]
CON_Score=[0,0,0,0,0,0,0,0,0,0]
#Konsensüs sekansı tutmak için boş dizi 
CON_SEKANS=[]

#matrix sütünlerini string diziye çevirme fonksyonu
def convert(s):
    new = ""
    for x in s:
        new += x  
    return new
#sütünler sayısı kader string dizi sekansı
seq1=convert(A[:,0])
seq2=convert(A[:,1])
seq3=convert(A[:,2])
seq4=convert(A[:,3])
seq5=convert(A[:,4])
seq6=convert(A[:,5])
seq7=convert(A[:,6])
seq8=convert(A[:,7])
seq9=convert(A[:,8])

          
# entropi hesabı
def entropy(dna_sequence):
    m = len(dna_sequence)
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
 
    value = 0
    for base in bases:
      
        n_i = bases[base]
        # n_i (# residues type i) / M (# residues in column)
        p_i = n_i / float(m)
        entropy_i = p_i * (math.log(p_i, 2))
        value += entropy_i
 
    return value * -1


#her sütün dizisinde elemanlar ikişer alıp ve distance mantığına göre karşılaşacaklar
def distance(a, b):
    if a == b:
        return 0 
    elif a == "-" or b == "-":
        return 1 
    else:
        return 2
#her sütün dizisinde elemanlar ikişer alıp ve eşititlik ve varlık mantığına göre karşılaşacaklar
def w_sum(a, b):
    
    if a == b:
        return 3 
    elif a == "-" or b == "-":
        return -1 
    else:
        return 2 
#her sütünün w_sum karşılaştımasını ve elde edilen sayısal sonuçların diziye aktarma fonksyonu    
def Som_Of_Pairs(seq,col):   
    n=len(seq)   
    sp=0
    
    for x in range(n):
        k=1
        
        for y in range(n-x-1):
            print (str(col)+".sutun sırasıyla edit uzaklık: "+ str(w_sum(seq[x], seq[x+k])))
            sp+=w_sum(seq[x], seq[x+k])
            k+=1   
    SOP_Score[col]=sp

def Sop_Score():
    score=0
    for x in range(len(SOP_Score)):
        score+=SOP_Score[x]
    return score    
        
        

def Consensus(seq,col):   
    n=len(seq)   
    sp=0
    
    for x in range(n):
        k=1
        
        for y in range(n-x-1):
            print ( str(col)+".sutun sırasıyla Konsensüs sekansları : "+ str(distance(seq[x], seq[x+k])))
            sp+=distance(seq[x], seq[x+k])
           
            k+=1

    CON_Score[col]=sp

def Con_Score():
    score=0
    for x in range(len(CON_Score)):
        score+=CON_Score[x]
    return score 

#her satır için en tekrarlı elemanı bulma       
def Con_Sekans_Value(seq):
    #varlık değerleri
    a=0
    c=0
    g=0
    t=0
    bos=0
    #aktarılan sonuç değer kümesi
    val=[''] 
    #varlık değişimi
    for x in range(len( seq)):
        if seq[x]=="A":
            a+=1
        elif seq[x] == "C":
            c+=1
        elif seq[x] == "G":
            g+=1
        elif seq[x] == "T":
            t+=1 
        else:
            bos+=1 
    #varlık sayısını tutmak için oluşturulan dizi        
    sek=[a,c,g,t,bos]
    #max index bulam
    x=sek.index(max(sek)) 
    #değer atama ve aktarma
    if x==0:
        val.append('A')
        return val[1] 
    elif x==1:
        val.append('C')
        return val[1] 
    elif x==2:
       val.append('G')
       return val[1] 
    elif x==3:
        val.append('T')
        return val[1] 
    
    else:
        val.append('-')
        return val[1] 
   
      
#sonuç Konsensüs sekansı nı oluşturan fonksyon
def Con_Sek_Dizi(seq):
        CON_SEKANS.append(Con_Sekans_Value(seq))
        
Con_Sek_Dizi(seq1)      
Con_Sek_Dizi(seq2)
Con_Sek_Dizi(seq3)
Con_Sek_Dizi(seq4)   
Con_Sek_Dizi(seq5)   
Con_Sek_Dizi(seq6)   
Con_Sek_Dizi(seq7)  
Con_Sek_Dizi(seq8)   
Con_Sek_Dizi(seq9) 

         
print("----------ENTROPY-------------")   
print("1.sütün dizisi= "+str(seq1)+"  Entropy=> "+str(entropy(seq1)))
print("2.sütün dizisi= "+str(seq2)+"  Entropy=> "+str(entropy(seq2)))
print("3.sütün dizisi= "+str(seq3)+"  Entropy=> "+str(entropy(seq3)))
print("4.sütün dizisi= "+str(seq4)+"  Entropy=> "+str(entropy(seq4)))
print("5.sütün dizisi= "+str(seq5)+"  Entropy=> "+str(entropy(seq5)))
print("6.sütün dizisi= "+str(seq6)+"  Entropy=> "+str(entropy(seq6)))
print("7.sütün dizisi= "+str(seq7)+"  Entropy=> "+str(entropy(seq7)))
print("8.sütün dizisi= "+str(seq8)+"  Entropy=> "+str(entropy(seq8)))
print("9.sütün dizisi= "+str(seq9)+"  Entropy=> "+str(entropy(seq9)))     
print("----------SUM OF PAIRS-------------")
print(Som_Of_Pairs(seq1,1))
print("1.sutun edit distance toplam: "+str(SOP_Score[1]))
print(Som_Of_Pairs(seq2,2))
print("2.sutun edit distance toplam: "+str(SOP_Score[2]))
print(Som_Of_Pairs(seq3,3))
print("3.sutun edit distance toplam: "+str(SOP_Score[3]))
print(Som_Of_Pairs(seq4,4))
print("4.sutun edit distance toplam: "+str(SOP_Score[4]))
print(Som_Of_Pairs(seq5,5))
print("5.sutun edit distance toplam: "+str(SOP_Score[5]))
print(Som_Of_Pairs(seq6,6))
print("6.sutun edit distance toplam: "+str(SOP_Score[6]))
print(Som_Of_Pairs(seq7,7))
print("7.sutun edit distance toplam: "+str(SOP_Score[7]))
print(Som_Of_Pairs(seq8,8))
print("8.sutun edit distance toplam: "+str(SOP_Score[8]))
print(Som_Of_Pairs(seq9,9))
print("9.sutun edit distance toplam: "+str(SOP_Score[9]))
print("--------------------------")
print ("toplam Sum Of Pairs scoru:  "+str(Sop_Score()))
print("----------CONSENSUS MAT-------------")
print(Consensus(seq1,1))
print("1.sutun Konsensüs toplam: "+str(CON_Score[1]))
print(Consensus(seq2,2))
print("2.sutun Konsensüs toplam: "+str(CON_Score[2]))
print(Consensus(seq3,3))
print("3.sutun Konsensüs toplam: "+str(CON_Score[3]))
print(Consensus(seq4,4))
print("4.sutun Konsensüs toplam: "+str(CON_Score[4]))
print(Consensus(seq5,5))
print("5.sutun Konsensüs toplam: "+str(CON_Score[5]))
print(Consensus(seq6,6))
print("6.sutun Konsensüs toplam: "+str(CON_Score[6]))
print(Consensus(seq7,7))
print("7.sutun Konsensüs toplam: "+str(CON_Score[7]))
print(Consensus(seq8,8))
print("8.sutun Konsensüs toplam: "+str(CON_Score[8]))
print(Consensus(seq9,9))
print("9.sutun Konsensüs toplam: "+str(CON_Score[9]))
print("--------------------------")
print ("toplam Konsensüs scoru:  "+str(Con_Score()))
print("----------CONSENSUS SERIES-------------")

print(CON_SEKANS)