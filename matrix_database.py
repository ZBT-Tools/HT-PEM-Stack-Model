
# coding: utf-8

# In[3]:

from scipy.linalg import block_diag
import numpy as np


# In[11]:

def element_to_node(n):#element to node, last row still need to be applied (.append())
    m = np.full((n,n),0.)
    for j in range (n):
        for i in range (n):
            if i is j and j >0 : 
                m[j,i] = 0.5
                m[j,i-1] = 0.5
    return m

def backward_matrix(n):#element backward
    m = np.full((n,n),0.)
    for j in range (n):
        for i in range (n):
            if i is n-1: m[j,i] = 0.5
            if i < n-1 and i>=j: m[j,i] = 1.
    return m

def forward_matrix(n):#element forward
    m = np.full((n,n),0.)
    for j in range (n):
        for i in range (n):
            if i is 0: m[j,i] = 0.5
            if i >= 1 and i<=j: m[j,i] = 1.
    return m

def node_to_element(n):
    m = np.full((n-1,n),0.)
    for j in range (n-1):
        for i in range (n):
            if j is i:
                m[j,i]= 0.5 
                m[j,i+1]= 0.5
    return m

def element_to_node_func(start,end,element_vector):
    node_vector = np.matmul(element_to_node(len(element_vector)),element_vector).tolist()
    node_vector[0] = start
    node_vector.append(end)
    return np.asarray(node_vector)

def temperature_matrix(n,M,mu_p,mu_g,mu_m,a):
    m = np.full((M*n*5,M*n*5),0.)
    w = 0
    while w<n*5:
        m[0+w,0+w] = -mu_g
        m[0+w,1+w] = mu_p+mu_g
        m[0+w,2+w] = -mu_p
        m[1+w,0+w] = mu_g+mu_m
        m[1+w,1+w] = -mu_g
        m[1+w,3+w] = -mu_m
        m[2+w,0+w] = -mu_m
        m[2+w,3+w] = mu_g + mu_m
        m[2+w,4+w] = -mu_g
        m[3+w,2+w] = -mu_p
        m[3+w,3+w] = -mu_g 
        m[3+w,4+w] = mu_p + mu_g
        m[4+w,1+w] = mu_p
        m[4+w,2+w] = -2.*mu_p -a
        m[4+w,4+w] = mu_p
        w =w+5
    return m

def b(n,M,d_x):
    m = np.full((n,n),0.)
    for j in range (n):
        for l in range (n):
            if l==0 and j==0:
                m[l,j] = -1.
                m[l,j+1] = 1.
            elif ((j==l and (j<n-1 and l <n)) and (j>0 and l>0)) and m[l,j]==0:
                m[l,j] = -2.
                m[l,j-1] = 1.
                m[l,j+1] = 1.
            elif j==n-1 and l==n-1:
                m[l,j] = 1.
                m[l,j-1] = -1.

    ####Building up the B ((n+1)*M)x((n+1)*M) sec. order ode Matrix based on zdm2         
    b = [m]*M
    B = block_diag(*b)/d_x**2.
    return B

def c(n,M,d_x):
    m = np.full((n*M,n*M),0.)
    for j in range (n*M):
        for l in range (n*M):
            if j==l  and (j>=n and j<M*n-n):
                m[l,j] = -2.
                m[l,j-(n)] = 1.
                m[l,j+(n)] = 1.
    return m/d_x**2.
