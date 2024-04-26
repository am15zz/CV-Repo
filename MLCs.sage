x, y, a = var('x y a')
half_field_size=25
equation1 = y==(1000/half_field_size)*x+1000
equation2 = (170)^2==(x+a+170)^2+(y-7.5-644.2)^2
p1,p2 = solve([equation1,equation2],x,y)
euclidean_distance=sqrt((p2[0].rhs()-p1[0].rhs())^2+(p2[1].rhs()-p1[1].rhs())^2)
MV6_mudrho = 0.0284334454335
MV10_mudrho = 0.0293966014917
MV18_mudrho = 0.032122200251
rho = 18
TVL=log(2)/(rho*(MV18_mudrho))
position=solve([euclidean_distance==TVL],a)
numeric_position=[value.rhs().numerical_approx() for value in position]
print(numeric_position)

 

