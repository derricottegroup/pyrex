import numpy as np
import matplotlib.pyplot as plt

N = 4
elst = (-5.21, -8.91, -9.78, -9.46)
exch = (13.14, 13.41, 14.01, 14.39)
indAB = (-1.20, -1.02, -1.19, -1.44)
indBA = (-6.36, -13.40, -12.93, -12.07)
disp = (-2.03, -2.05, -2.77, -3.07)
total = (-1.65, -11.97, -12.66, -11.65)
men_means = (20, 35, 30, 35, 27)
women_means = (25, 32, 34, 20, 25)

ind = np.arange(N) 
width = 0.1       
plt.bar(ind, elst, width, label='elst', color="red")
plt.bar(ind + width, exch, width,
    label='exch', color="green")
plt.bar(ind + 2.0*width, indAB, width,
    label='indAB', color="blue")
plt.bar(ind + 3.0*width, indBA, width,
    label='indBA', color="orange")
plt.bar(ind + 4.0*width, disp, width,
    label='disp', color="magenta")
plt.bar(ind + 5.0*width, total, width,
    label='total', color="black")


plt.axhline(y=0,color='black',linestyle="--",linewidth=2.5)

plt.ylabel('$w_1 + w_2$ (kcal mol$^{-1}$')
plt.title('CX$_3$ --- OH')

plt.xticks(ind + 5.0*width / 2, ('X=H', 'X=F', 'X=Cl', 'X=Br'))
plt.legend(loc='best')
plt.show()
