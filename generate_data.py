import random
a = open('non_over_lap_area','w')
for i in range(0,333):
	a.write(str(random.uniform(2,6))+"  "+str(random.uniform(0,4))+"\n")
	
for i in range(0,333):
	a.write(str(random.uniform(10,12))+"  "+str(random.uniform(20,22))+"\n")

for i in range(0,334):
	a.write(str(random.uniform(14,18))+"  "+str(random.uniform(30,35))+"\n")
