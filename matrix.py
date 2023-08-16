import numpy as np
import pandas as pd
from sympy import Matrix
import subprocess
import matplotlib.pyplot as plt
import itertools

def to_mathematica(matrix):
    if isinstance(matrix, np.ndarray):
        matrix = matrix.tolist()
    
    rows = ['{' + ', '.join(map(str, row)) + '}' for row in matrix]
    return '{' + ', '.join(rows) + '}'

def read_output(energy_array, print_output=False, print_mathematica=False):
    np.set_printoptions(precision=3)

    output_data = pd.read_csv("./output.txt",delimiter='\t',header=1)
    functions_numbered = []
    nmin = abs(output_data.iloc[0:,0].min())
    nmax = abs(output_data.iloc[0:,0].max())

    umin = abs(output_data.iloc[0:,2].min())
    umax = abs(output_data.iloc[0:,2].max())

    mmin = abs(output_data.iloc[0:,4].min())
    mmax = abs(output_data.iloc[0:,4].max())

    for i in range(nmin,nmax+1):
        for j in range(umin,umax+1):
            for k in range(mmin,mmax+1):
                functions_numbered.append(np.array([i,j,k]))

    functions_numbered = np.array(functions_numbered)

    matrix = np.zeros((len(functions_numbered),len(functions_numbered)))

    for row in output_data.iloc: #daje obiekty w postaci rzędów
        function_1 = np.array([row.n1,row.u1,row.m1])
        function_2 = np.array([row.n2,row.u2,row.m2])

        index_1 = np.where(np.all(functions_numbered==function_1,axis=1))[0]
        index_2 = np.where(np.all(functions_numbered==function_2,axis=1))[0]
        
        matrix[index_2, index_1] += row.integral_value

    m = Matrix(matrix)
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    diagonal_matrix = np.diag(eigenvalues)

    energy_array.append(np.min(eigenvalues))
    """
    if not np.allclose(matrix, np.transpose(matrix),1e-14):
        print(matrix-np.transpose(matrix))
        raise Exception("energy matrix doesn't seem to be symmetric")
    """
    if print_output:
        print("order of functions (n_1 = -n_2) by multiindex:\n [n,u,m]")
        print(functions_numbered)
        print("energy matrix:")
        print(matrix)
        print("diagonalized energy matrix:")
        print(diagonal_matrix)
    if print_mathematica:
        print(to_mathematica(matrix))
    #print(f"kolejność liczb kwantowych macierzy (n,l,m): \n\n{functions_numbered}\n\n macierz: \n\n{matrix}\n\n macierz po diagonalizacji: \n\n{diagonal_matrix}")
    return

def read_diagonal_output(energy_array):
    np.set_printoptions(precision=3)

    output_data = pd.read_csv("./output.txt",delimiter='\t',header=1)
    functions_numbered = []
    matrix = np.zeros((len(output_data),len(output_data)))
    for i in enumerate(output_data.iloc[:,6]):
        matrix[i[0]][i[0]] += i[1]
    m = Matrix(matrix)
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    diagonal_matrix = np.diag(eigenvalues)

    energy_array.append(np.min(eigenvalues))
    #print(f"kolejność liczb kwantowych macierzy (n,l,m): \n\n{functions_numbered}\n\n macierz: \n\n{matrix}\n\n macierz po diagonalizacji: \n\n{diagonal_matrix}")
    return

def generate_output_for_given_g(g, l, nmin, nmax, mmin, mmax, umin, umax):
    subprocess.run(["./g-plot", "--change_g", f"{g}", "--change_l", f"{l}","--minmax", f"{nmin}", f"{nmax}", f"{mmin}", f"{mmax}", f"{umin}", f"{umax}"])
    return

def plot_array(x, y):
    x = np.array(x)
    y = np.array(y)
    plt.plot(x,y,'o', ls='-')
    plt.show()
    return

def generate_single_matrix(g,length,nmin, nmax, mmin, mmax, umin, umax):
    energy_array = []
    generate_output_for_given_g(g, length, nmin, nmax, mmin, mmax, umin, umax)
    read_output(energy_array,True)


def generate_many_minimum_energy_points(nmin, nmax, mmin, mmax, umin, umax):
    energy_array = [] # lowest energy of a system
    x_array = []

    #l - 0->10
    #g - 0->2000
    max_changing_value = 1
    for x in np.arange(0.001,max_changing_value,max_changing_value/20):

        generate_output_for_given_g(1, x, nmin, nmax, mmin, mmax, umin, umax)
        x_array.append(x)
        read_output(energy_array)

    print(energy_array)
    plot_array(x_array, energy_array)

def find_most_probable_minmax(length, num, nmin=0, umin=0, mmin=0):#num is number of matrix rows and columns
    all_values = []
    all_indices = []
    def energy_func(n,u,m):
        return n**2*np.pi**2*4/length**2 + 2*u + m
    # Generate all possible combinations for indices (you can adjust the range based on your needs)
    for i, j, k in itertools.product(range(10), repeat=3): # Here I use a range of 10 for each index as an example.
        value = energy_func(i, j, k)
        all_values.append(value)
        all_indices.append([i,j,k])

    sorted_data = sorted(zip(all_values, all_indices))
    sorted_vals, sorted_indices = zip(*sorted_data)
    #print(sorted_indices[:num])
    srtd=sorted_indices[:num]
    max_values = [max(values) for values in zip(*srtd)]
    #print(max_values)
    nmax,umax,mmax =max_values
    return nmax,umax,mmax
    
    


def main():
    g=1000
    length=1
    nmax, umax, mmax = find_most_probable_minmax(length, 5)
    print("nmax, umax, mmax: ",find_most_probable_minmax(length, 5))
    nmin, umin, mmin = 0,0,0
    #nmin, nmax, mmin, mmax, umin, umax = 0, 0, 0, 0, 0, 0
    generate_single_matrix(g,length,nmin, nmax, mmin, mmax, umin, umax)
    return 0

if __name__ == "__main__":
    main()
    
    #narazie się nie zgadza, jednak może to wynikać bezpośrednio stąd, że mamy dwie różne liczby kwantowe - michał napisał kod dla dwóch PRZECIWNYCH
    #żeby rzecz przyspieszyć po prostu zachowam macierz energii delty i przemnożę je przez g
    #przyjżę się też sposobowi w jaki powinienem liczyć macierz energii potencjalnej dla drugiej (i wyższej) kwantyzacji, ale to wymaga dłuższej lektury 

    
