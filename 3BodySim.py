#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: JasonQSY <jasonsyqian@gmail.com>
# License: MIT

import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt, sin, cos

# d1 jupyter
d1 = (1.898, 0, -6.286, -13.6 * 7, 0)
# d2 satellite
#d2 = (0.0005, 3.082, 0, 0, 22.5 * 7)
d2 = (0.0005, 1.482, 0, 0, 22.5 * 12)

d3=(4e-4,10)
d4=(2000.0)

clocktick=0 # 0 for fastest animation, microseconds elapsed between steps

m1,x01,y01,vx01,vy01 = d1
m2,x02,y02,vx02,vy02 = d2
h,tf = d3
M = d4

G=-4*pi**2

# E = K1 + V1 + K2 + V2 + V12
E0 = m1*(0.5*(vx01**2 + vy01**2) + G*M/sqrt(x01**2 + y01**2))
E0+= m2*(0.5*(vx02**2 + vy02**2) + G*M/sqrt(x02**2 + y02**2))
E0+= G*m1*m2/sqrt((x01-x02)**2+(y01-y02)**2)

x1= [x01]; y1= [y01]; vx1=[vx01]; vy1=[vy01]
x2= [x02]; y2= [y02]; vx2=[vx02]; vy2=[vy02]
E=[E0]; e=[0.0]

i=0
ax1=[ G*M*x1[i]/(x1[i]**2+y1[i]**2)**1.5 + G*m2*(x1[i]-x2[i])/((x1[i]-x2[i])**2+(y1[i]-y2[i])**2)**1.5]
ay1=[ G*M*y1[i]/(x1[i]**2+y1[i]**2)**1.5 + G*m2*(y1[i]-y2[i])/((x1[i]-x2[i])**2+(y1[i]-y2[i])**2)**1.5]
ax2=[ G*M*x2[i]/(x2[i]**2+y2[i]**2)**1.5 + G*m1*(x2[i]-x1[i])/((x2[i]-x1[i])**2+(y2[i]-y1[i])**2)**1.5]
ay2=[ G*M*y2[i]/(x2[i]**2+y2[i]**2)**1.5 + G*m1*(y2[i]-y1[i])/((x2[i]-x1[i])**2+(y2[i]-y1[i])**2)**1.5]

t=np.arange(0,tf,h)
for i in xrange(0,len(t)-1):
	x1.append(  x1[i]  + h*vx1[i] + h*h*ax1[i]/2     )
	y1.append(  y1[i]  + h*vy1[i] + h*h*ay1[i]/2     )
	x2.append(  x2[i]  + h*vx2[i] + h*h*ax2[i]/2     )
	y2.append(  y2[i]  + h*vy2[i] + h*h*ay2[i]/2     )

	ax1.append( G*M*x1[i+1]/(x1[i+1]**2+y1[i+1]**2)**1.5 +  G*m2*(x1[i+1]-x2[i+1])/((x1[i+1]-x2[i+1])**2+(y1[i+1]-y2[i+1])**2)**1.5 )
	ay1.append( G*M*y1[i+1]/(x1[i+1]**2+y1[i+1]**2)**1.5 +  G*m2*(y1[i+1]-y2[i+1])/((x1[i+1]-x2[i+1])**2+(y1[i+1]-y2[i+1])**2)**1.5 )
	ax2.append( G*M*x2[i+1]/(x2[i+1]**2+y2[i+1]**2)**1.5 +  G*m1*(x2[i+1]-x1[i+1])/((x2[i+1]-x1[i+1])**2+(y2[i+1]-y1[i+1])**2)**1.5 )
	ay2.append( G*M*y2[i+1]/(x2[i+1]**2+y2[i+1]**2)**1.5 +  G*m1*(y2[i+1]-y1[i+1])/((x2[i+1]-x1[i+1])**2+(y2[i+1]-y1[i+1])**2)**1.5 )

	vx1.append( vx1[i] + h*(ax1[i+1] + ax1[i] )/2 )
	vy1.append( vy1[i] + h*(ay1[i+1] + ay1[i] )/2 )
	vx2.append( vx2[i] + h*(ax2[i+1] + ax2[i] )/2 )
	vy2.append( vy2[i] + h*(ay2[i+1] + ay2[i] )/2 )

	Etmp = m1*(0.5*(vx1[i]**2 + vy1[i]**2) + G*M/sqrt(x1[i]** 2+ y1[i]**2))
	Etmp+= m2*(0.5*(vx2[i]**2 + vy2[i]**2) + G*M/sqrt(x2[i]** 2+ y2[i]**2))
	Etmp+= G*m1*m2/sqrt((x1[i]-x2[i])**2+(y1[i]-y2[i])**2)

	E.append( Etmp )
	e.append( 100*( Etmp - E0 )/E0 )

plt.figure(0)
plt.plot(x1, y1, 'r')

x3 = [x2[0]]
y3 = [y2[0]]
for i in range(1, 100):
	x3.append(x3[0] / i)
	y3.append(y3[0] / i)

plt.plot(x3, y3, 'b')
plt.plot(0,0,'ro', 1.496, 0, 'go')
plt.xlim(-7, 7)
plt.ylim(-7, 4)
plt.grid(True)
plt.xlabel(r'$x(10^9 km)$')
plt.ylabel(r'$y(10^9 km)$')
plt.legend(('Jupiter', 'waste', 'sun', 'earth'))
plt.show()
