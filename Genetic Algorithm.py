
# coding: utf-8

# In[14]:
import random
import math
import matplotlib.pyplot as plt 


# In[15]:

# population size
popSize = 500
# chromosome size
chromSize = 40
# max value of a chromosome
maxValue = math.pi
speed = 100
gravity = 9.8
generation = 2000
matingPos = 0.6
mutationPos = 0.1


# In[16]:

def GenerateChromosome(popSize, chromSize):
    pop = []
    for i in range(popSize):
        pop.append([])
        for j in range(chromSize):
            pop[i].append(random.randint(0,1))
    return pop


# In[17]:

def DecodeChromsome(pop,chromSize,maxValue,speed):
    radian = []
    decodeVal = []        
    for i in pop:
        r = ChromsomeToRadian(i,chromSize,maxValue)
        distance = CalculateDistance(r, speed)
        radian.append(r)
        decodeVal.append(distance)
    return radian,decodeVal

def CalculateDistance(rad,speed):
    ver = speed * math.sin(rad)
    hor = speed * math.cos(rad)
    time = (ver / gravity) * 2
    dis = hor * time
    return dis
    


# In[18]:

# 消除負值
# 不將負值設為0的原因是因為若所有基因皆為0, 分母將除以0
def Eliminate(value):
    for i in range(len(value)):
        if(value[i] < 1):
            value[i] = 1 


# In[19]:

def ChromsomeToRadian(gene,chromSize,maxValue):
    temp = 0
    chromMaxValue = pow(2,chromSize) - 1
    
    for i in range(len(gene)):
        temp = temp*2 + gene[i]
    temp = temp / chromMaxValue * maxValue
    return temp

def RadianToDegree(radian):
    return radian * 180 / math.pi


# In[20]:

def Best(pop,radian,distance):
    dis = 0
    degree = 0    
    gene = []
    for i in range(len(radian)):
        if(dis < distance[i]):
            dis = distance[i]
            degree = RadianToDegree(radian[i])            
            gene = pop[i][:]
    return gene,degree,dis


# In[21]:

# 計算累積機率
def CumSum(pos):
    temp = 0
    for i in range(len(pos) - 1):
        temp += pos[i]
        pos[i] = temp
    pos[len(pos)-1] = 1
    

# 表現越好的基因佔有比例越大
def Selection(pop,val):    
    totalVal = sum(val)
    pos = []
    for i in val:
        pos.append(i/totalVal)           
    CumSum(pos)     
    i = 0
    idx = 0
    newPop = [[]]
    while(i < len(pop)):              
        while(i < round(pos[idx] * len(pop))):            
            newPop.append(pop[idx][:])            
            i +=1                        
        idx+=1                        
    return newPop[1:]    


# In[22]:

# 交配
def Mating(pop):    
    for i in range(len(pop)):
        if(random.random() <= matingPos):            
            idx = random.randint(0,len(pop) - 2)
            if idx >= i:
                idx+=1
            mPoint = random.sample(range(0,len(pop[0])),2)
            mPoint.sort()
            pop[i][mPoint[0]:mPoint[1]] = pop[idx][mPoint[0]:mPoint[1]]     

# 突變
def Mutation(pop):
    for i in range(len(pop)):        
        if(random.random() <= mutationPos):  
            mPoint = random.randint(0,len(pop[0])-1)            
            if(pop[i][mPoint] == 0):
                pop[i][mPoint] = 1                        
            else:
                pop[i][mPoint] = 0                                            
           


# In[23]:

pop = GenerateChromosome(popSize, chromSize)
results = [[]]

bestDistance = 0
bestDegree = 0
bestGene = []
for i in range(generation):        
    radian,val = DecodeChromsome(pop,chromSize,maxValue,speed)
    Eliminate(val)    
    bestGene,degree, distance = Best(pop,radian,val)
    if(distance > bestDistance):
        bestDistance = distance
        bestDegree = degree        
    results.append([degree,distance])
    pop = Selection(pop,val)        
    Mating(pop)            
    Mutation(pop) 
    

results = results[1:]
print("Best Degree : " + str(bestDegree) + " Best Distance : " + str(bestDistance) + " Best Gene : " + str(bestGene))

# Plot
X = []
Ydegree = []
Ydistance = []
for i in range(generation):
    X.append(i)        
    Ydegree.append(results[i][0])
    Ydistance.append(results[i][1])



# In[25]:

plt.rcParams['figure.figsize'] = [20, 15]

plt.figure()
plt.plot(X, Ydegree)
plt.show()

plt.rcParams['figure.figsize'] = [20, 15]

plt.figure()
plt.plot(X, Ydistance)
plt.show()


