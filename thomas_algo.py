

import numpy as np
import matplotlib.pyplot as plt

def A(x): return (1.0 / x)
def B(x): return (0.0)
def C(x): return (1 / x **2)

def thomas(a, b, c, d):
    n = len(d)
    c1 = [0]*n
    d1 = [0]*n
    y = [0]*n

    c1[0] = c[0] / b[0]
    d1[0] = d[0] / b[0]

    for i in range(1, n, 1):
        c1[i] = c[i] / (b[i] - a[i] * c1[i - 1])
        d1[i] = (d[i] - a[i] * d1[i - 1]) / (b[i] - a[i] * c1[i - 1])

    y[n - 1] = d1[n - 1]
    for i in range(n - 2, -1, -1):
        y[i] = d1[i] - c1[i] * y[i + 1]

    return y


def helper(a_, b_, h, y_a, y_b):
    n = int((b_ - a_) / h)+1
    a = [0]*(n-1)
    b = [0]*(n-1)
    c = [0]*(n-1)
    d = [0]*(n-1)
    for i in range(1, n):
        x = a_ + i * h
        a[i-1] = (1.0 / (h ** 2))- (A(x) / (2.0 * h))
        b[i-1] = (-2.0 / (h ** 2)) + B(x)
        c[i-1] = (1.0 / (h ** 2)) + (A(x) / (2.0 * h))
        if i == 1:
            d[i-1] = C(x) - a[i-1] * y_a
        elif i == n - 1:
            d[i-1] = C(x) - c[i-1] * y_b
        else:
            d[i-1] = C(x)

    return [y_a] + thomas(a, b, c, d) + [y_b]


def Output(x , y, h):
    print("\t\tOutput for h="+str(h)+"\n\n")
    print("\t\tX\t\t\tY\n\n")
    for i in range(len(x)):
        print("\t\t"+str('%.3f'%x[i])+"\t\t"+str(y[i])+"\n")

def main():
    a = 1
    b = 1.4
    h = [0.1, 0.05, 0.001]
    y_a = 0
    y_b = 0.0566
    n1=int((b - a) / h[0])+1
    n2=int((b - a) / h[1])+1
    n3=int((b - a) / h[2])+1
    y_1 = [0]*(n1+1)
    y_2 = [0]*(n2+1)
    y_3 = [0]*(n3+1)
    x_1= np.linspace(a, b, int((b - a) / h[0]+1) + 1)
    x_2= np.linspace(a, b, int((b - a) / h[1]+1) + 1)
    x_3= np.linspace(a, b, int((b - a) / h[2]+1) + 1)
    y_1 = helper(a, b, h[0], y_a, y_b)
    y_2 = helper(a, b, h[1], y_a, y_b)
    y_3 = helper(a, b, h[2], y_a, y_b)
        
    Output(x_1,y_1,h[0])
    Output(x_2,y_2,h[1])
    Output(x_3,y_3,h[2])

    p1,p2,p3=plt.plot(x_3,np.interp(x_3,x_1,y_1),'r',x_3,np.interp(x_3,x_2,y_2),'g',x_3,y_3,'b')
    plt.legend([p1, p2, p3], ["h = 0.1", "h =0.05", "h = 0.001"], loc =4)
    plt.ylabel("Y(x)")
    plt.xlabel ("X")
    plt.show()

main()

import numpy as np
import matplotlib.pyplot as plt

#Given Differential Equation: y" = x + y   y(0) = y(1) = 0

def A(x): return (0.0)
def B(x): return (-1)
def C(x): return (x)

def thomas(a, b, c, d):
    n = len(d) -1
    c1 = [0]*n
    d1 = [0]*n
    y = [0]*n

    c1[0] = c[0] / b[0]
    d1[0] = d[0] / b[0]

    for i in range(1, n, 1):
        c1[i] = c[i] / (b[i] - a[i] * c1[i - 1])
        d1[i] = (d[i] - a[i] * d1[i - 1]) / (b[i] - a[i] * c1[i - 1])

    y[n - 1] = d1[n - 1]
    for i in range(n - 2, -1, -1):
        y[i] = d1[i] - c1[i] * y[i + 1]

    return y


def helper(a_, b_, h, y_a, y_b):
    n = int((b_ - a_) / h)+1
    a = [0]*(n-1)
    b = [0]*(n-1)
    c = [0]*(n-1)
    d = [0]*(n-1)
    for i in range(1, n):
        x = a_ + i * h
        a[i-1] = (1.0 / (h ** 2))- (A(x) / (2.0 * h))
        b[i-1] = (-2.0 / (h ** 2)) + B(x)
        c[i-1] = (1.0 / (h ** 2)) + (A(x) / (2.0 * h))
        if i == 1:
            d[i-1] = C(x) - a[i-1] * y_a
        elif i == n - 1:
            d[i-1] = C(x) - c[i-1] * y_b
        else:
            d[i-1] = C(x)

    return [y_a] + thomas(a, b, c, d) + [y_b]


def Output(x , y, h):
    print("\t\tOutput for h="+str(h)+"\n\n")
    print("\t\tX\t\t\tY\n\n")
    for i in range(len(x)):
        print("\t\t"+str('%.3f'%x[i])+"\t\t"+str(y[i])+"\n")

def main():
    a = 0
    b = 1.0
    h = [0.1, 0.05, 0.001]
    y_a = 0
    y_b = 0
    n1 =int((b - a) / h[0])+1
    n2=int((b - a) / h[1])+1
    n3=int((b - a) / h[2])+1
    y_1 = [0]*(n1)
    y_2 = [0]*(n2)
    y_3 = [0]*(n3)
    x_1= np.linspace(a, b, int((b - a) / h[0]) + 1)
    x_2= np.linspace(a, b, int((b - a) / h[1]) + 1)
    x_3= np.linspace(a, b, int((b - a) / h[2]) + 1)
    y_1 = helper(a, b, h[0], y_a, y_b)
    y_2 = helper(a, b, h[1], y_a, y_b)
    y_3 = helper(a, b, h[2], y_a, y_b)
        
    Output(x_1,y_1,h[0])
    Output(x_2,y_2,h[1])
    Output(x_3,y_3,h[2])

    p1,p2,p3=plt.plot(x_3,np.interp(x_3,x_1,y_1),'r',x_3,np.interp(x_3,x_2,y_2),'g',x_3,y_3,'b')
    plt.legend([p1, p2, p3], ["h = 0.1", "h =0.05", "h = 0.001"], loc =4)
    plt.ylabel("Y(x)")
    plt.xlabel ("X")
    plt.show()

main()

import numpy as np
import matplotlib.pyplot as plt

# Given Differential Equation : y"-2y = 0 ,    y(0) = 1, & y'(1) = 0
# Use 2nd order backward difference for y'(1) = 0 

def A(x): return 0
def B(x): return -2.0
def C(x): return 0

def thomas(a, b, c, d):
    n = len(d)
    c1 = [0]*n
    d1 = [0]*n
    y = [0]*n

    c1[0] = c[0] / b[0]
    d1[0] = d[0] / b[0]

    for i in range(1, n, 1):
        c1[i] = c[i] / (b[i] - a[i] * c1[i - 1])
        d1[i] = (d[i] - a[i] * d1[i - 1]) / (b[i] - a[i] * c1[i - 1])

    y[n - 1] = d1[n - 1]
    for i in range(n - 2, -1, -1):
        y[i] = d1[i] - c1[i] * y[i + 1]

    return y


def helper(a1, b1, c1, h, a_, b_):
    n = int((b_ - a_) / (1.0*h))
    a = [0]*(n-1)
    b = [0]*(n-1)
    c = [0]*(n-1)
    d = [0]*(n-1)
    for i in range(1, n):
        x = a_ + i * h
        a[i-1] = (1.0 / (h ** 2))- (A(x) / (2.0 * h))
        b[i-1] = (-2.0 / (h ** 2)) + B(x)
        c[i-1] = (1.0 / (h ** 2)) + (A(x) / (2.0 * h))
        d[i-1] = C(x)
        if i == 1:
            den0 = a1[0] - (1.5 * b1[0] / h);
            b[i-1] += a[i-1] * (-2 * b1[0] / h) / den0;
            c[i-1] += a[i-1] * (0.5 * b1[0] /h) / den0;
            d[i-1] -= a[i-1] * (c1[0] / den0);


    den1 = a1[1] + (1.5 * b1[1] / h);
    b[i-1] += c[i-1] * (2 * b1[1] / h) / den1;
    a[i-1] += c[i-1] * (-0.5 * b1[1] /h) / den1;
    d[i-1] -= c[i-1] * (c1[1]/ den1);
    y=thomas(a, b, c, d)
    y_a=((-2 * b1[0]/ h)*y[0] + (0.5 * b1[0] /h)*y[1] +c1[0])/den0;
    y_b= ((2 * b1[1] / h)*y[-1] +(-0.5 * b1[1] /h)*y[-2] + c1[1])/den1;

    return [y_a] + y + [y_b]


def Output(x , y, h):
    print("\t\tOutput for h="+str(h)+"\n\n")
    print("\t\tX\t\t\tY\n\n")
    for i in range(len(x)):
        print("\t\t"+str('%.3f'%x[i])+"\t\t"+str(y[i])+"\n")


def main():
    a = 0
    b = 1
    h = [0.1, 0.05, 0.001]
    a1= [1.0, 0]
    b1=[0, 1.0]
    c1=[1.0, 0]
    n1 = int((b - a) / h[0])+1
    n2 = int((b - a) / h[1])+1
    n3 = int((b - a) / h[2])+1
    y_1 = [0]*(n1+1)
    y_2 = [0]*(n2+1)
    y_3 = [0]*(n3+1)
    x_1 = np.linspace(a, b, int((b - a) / (1.0*h[0])+1))
    x_2 = np.linspace(a, b, int((b - a) / (1.0*h[1])+1))
    x_3 = np.linspace(a, b, int((b - a) / (1.0*h[2])+1))
    y_1 = helper(a1, b1, c1, h[0], a, b)
    y_2 = helper(a1, b1, c1, h[1], a, b)
    y_3 = helper(a1, b1, c1, h[2], a, b)

    Output(x_1,y_1,h[0])
    Output(x_2,y_2,h[1])
    Output(x_3,y_3,h[2])

    
    p1,p2,p3=plt.plot(x_3,np.interp(x_3,x_1,y_1),'r',x_3,np.interp(x_3,x_2,y_2),'g',x_3,y_3,'b')
    plt.legend([p1, p2, p3], ["h = 0.1", "h =0.05", "h = 0.001"], loc = 1)
    plt.ylabel("Y")
    plt.xlabel ("X")
    plt.show()

main()

import numpy as np
import matplotlib.pyplot as plt

# Given Differential Equation : y"-2xy'-2y = -4x ,    y(0) - y'(0) = 0, & 2y(1) - y'(1) = 1 

def A(x): return (-2.0*x)
def B(x): return (-2.0)
def C(x): return (-4.0*x)

def thomas(a, b, c, d):
    n = len(d)
    c1 = [0]*n
    d1 = [0]*n
    y = [0]*n

    c1[0] = c[0] / b[0]
    d1[0] = d[0] / b[0]

    for i in range(1, n, 1):
        c1[i] = c[i] / (b[i] - a[i] * c1[i - 1])
        d1[i] = (d[i] - a[i] * d1[i - 1]) / (b[i] - a[i] * c1[i - 1])

    y[n - 1] = d1[n - 1]
    for i in range(n - 2, -1, -1):
        y[i] = d1[i] - c1[i] * y[i + 1]

    return y


def helper(a1, b1, c1, h, a_, b_):
    n = int((b_ - a_) / (1.0*h))
    a = [0]*(n-1)
    b = [0]*(n-1)
    c = [0]*(n-1)
    d = [0]*(n-1)
    for i in range(1, n):
        x = a_ + i * h
        a[i-1] = (1.0 / (h ** 2))- (A(x) / (2.0 * h))
        b[i-1] = (-2.0 / (h ** 2)) + B(x)
        c[i-1] = (1.0 / (h ** 2)) + (A(x) / (2.0 * h))
        d[i-1] = C(x)
        if i == 1:
            den0 = a1[0] - (1.5 * b1[0] / h);
            b[i-1] += a[i-1] * (-2 * b1[0] / h) / den0;
            c[i-1] += a[i-1] * (0.5 * b1[0] /h) / den0;
            d[i-1] -= a[i-1] * (c1[0] / den0);


    den1 = a1[1] + (1.5 * b1[1] / h);
    b[i-1] += c[i-1] * (2 * b1[1] / h) / den1;
    a[i-1] += c[i-1] * (-0.5 * b1[1] /h) / den1;
    d[i-1] -= c[i-1] * (c1[1]/ den1);
    y=thomas(a, b, c, d)
    y_a=((-2 * b1[0]/ h)*y[0] + (0.5 * b1[0] /h)*y[1] +c1[0])/den0;
    y_b= ((2 * b1[1] / h)*y[-1] +(-0.5 * b1[1] /h)*y[-2] + c1[1])/den1;

    return [y_a] + y + [y_b]


def Output(x , y, h):
    print("\t\tOutput for h="+str(h)+"\n\n")
    print("\t\tX\t\t\tY\n\n")
    for i in range(len(x)):
        print("\t\t"+str('%.3f'%x[i])+"\t\t"+str(y[i])+"\n")


def main():
    a = 0
    b = 1
    h = [0.1, 0.05, 0.001]
    a1= [1.0,2.0]
    b1=[-1.0,-1.0]
    c1=[0.0,1.0]
    n1 = int((b- a) / h[0])+1
    n2 = int((b - a) / h[1])+1
    n3 = int((b - a) / h[2])+1
    y_1 = [0]*(n1+1)
    y_2 = [0]*(n2+1)
    y_3 = [0]*(n3+1)
    x_1 = np.linspace(a, b, int((b - a) / (1.0*h[0])+1))
    x_2 = np.linspace(a, b, int((b - a) / (1.0*h[1])+1))
    x_3 = np.linspace(a, b, int((b - a) / (1.0*h[2])+1))
    y_1 = helper(a1, b1, c1, h[0], a, b)
    y_2 = helper(a1, b1, c1, h[1], a, b)
    y_3 = helper(a1, b1, c1, h[2], a, b)

    Output(x_1,y_1,h[0])
    Output(x_2,y_2,h[1])
    Output(x_3,y_3,h[2])

    
    p1,p2,p3=plt.plot(x_3,np.interp(x_3,x_1,y_1),'r',x_3,np.interp(x_3,x_2,y_2),'g',x_3,y_3,'b')
    plt.legend([p1, p2, p3], ["h = 0.1", "h =0.05", "h = 0.001"], loc = 4)
    plt.ylabel("Y")
    plt.xlabel ("X")
    plt.show()

main()

import numpy as np
import matplotlib.pyplot as plt

#Given Differential Equation: y"+2xy'+2y = 4x ,   y(0) = 1 and y(0.5) = 1.279 

def A(x): return (2.0*x)
def B(x): return (2.0)
def C(x): return (4.0*x)

def thomas(a, b, c, d):
    n = len(d)
    c1 = [0]*n
    d1 = [0]*n
    y = [0]*n

    c1[0] = c[0] / b[0]
    d1[0] = d[0] / b[0]

    for i in range(1, n):
        c1[i] = c[i] / (b[i] - a[i] * c1[i - 1])
        d1[i] = (d[i] - a[i] * d1[i - 1]) / (b[i] - a[i] * c1[i - 1])

    y[n - 1] = d1[n - 1]
    for i in range(n - 2, -1, -1):
        y[i] = d1[i] - c1[i] * y[i + 1]

    return y


def helper(a_, b_, h, y_a, y_b):
    n = int((b_ - a_) / h)
    a = [0]*(n-1)
    b = [0]*(n-1)
    c = [0]*(n-1)
    d = [0]*(n-1)
    for i in range(1, n):
        x = a_ + i * h
        a[i-1] = (1.0 / (h ** 2)) - (A(x) / (2.0 * h))
        b[i-1] = (-2.0 / (h ** 2)) + B(x)
        c[i-1] = (1.0 / (h ** 2)) + (A(x) / (2.0 * h))
        if i == 1:
            d[i-1] = C(x) - a[i-1] * y_a
        elif i == n - 1:
            d[i-1] = C(x) - c[i-1] * y_b
        else:
            d[i-1] = C(x)

    return [y_a] + thomas(a, b, c, d) + [y_b]


def Output(x , y, h):
    print("\t\tOutput for h="+str(h)+"\n\n")
    print("\t\tX\t\t\tY\n\n")
    for i in range(len(x)):
        print("\t\t"+str('%.3f'%x[i])+"\t\t"+str(y[i])+"\n")

def main():
    a = 0
    b = 0.5
    h = [0.1, 0.05, 0.001]
    y_a = 1
    y_b = 1.279
    n1 = int((b - a) / h[0])+1
    n2 = int((b - a) / h[1])+1
    n3 = int((b - a) / h[2])+1
    y_1 = [0]*(n1)
    y_2 = [0]*(n2)
    y_3 = [0]*(n3)
    x_1 = np.linspace(a, b, int((b - a) / h[0])+1)
    x_2 = np.linspace(a, b, int((b - a) / h[1])+1)
    x_3 = np.linspace(a, b, int((b - a) / h[2])+1)
    y_1 = helper(a, b, h[0], y_a, y_b)
    y_2 = helper(a, b, h[1], y_a, y_b)
    y_3 = helper(a, b, h[2], y_a, y_b)

    Output(x_1,y_1,h[0])
    Output(x_2,y_2,h[1])
    Output(x_3,y_3,h[2])

    p1,p2,p3=plt.plot(x_3,np.interp(x_3,x_1,y_1),'r',x_3,np.interp(x_3,x_2,y_2),'g',x_3,y_3,'b')
    plt.legend([p1, p2, p3], ["h = 0.1", "h =0.05", "h = 0.001"], loc =4)
    plt.ylabel("Y(x)")
    plt.xlabel ("X")
    plt.show()

main()