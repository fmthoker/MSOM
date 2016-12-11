import random
a = open('non_over_lap_area','w')
for i in range(0,333):
	a.write(str(random.uniform(2,4))+"  "+str(random.uniform(0,4))+"\n")
	
for i in range(0,333):
	a.write(str(random.uniform(6,8))+"  "+str(random.uniform(8,10))+"\n")

for i in range(0,334):
	a.write(str(random.uniform(7,9))+"  "+str(random.uniform(0,3))+"\n")
