# min(2*x1 − 3*x2 −    5*x4)
#      −x1 +   x2 − x3 − x4 ≤ 8
#     2*x1 + 4*x2           ≤ 10
#                   x3 + x4 ≤ 3
#                        xi ≥ 0

def simplex_method(matrix):

    print(matrix)

    while(1):
        minColIndex = matrix[0].index(min(matrix[0][1 : ]))

        if(matrix[0][minColIndex] >= 0):
            break

        row = []
        for i in range(1, len(matrix)):
            if(matrix[i][0] != 0 and matrix[i][minColIndex] > 0):
                #print(matrix[i][0], " ", matrix[i][minColIndex], " ", matrix[i][0] / matrix[i][minColIndex])
                row.append(matrix[i][0] / matrix[i][minColIndex])
            else:
                row.append(float('inf'))

        minRowIndex = row.index(min(row)) + 1

        print(matrix[0][minColIndex], " ", matrix[minRowIndex][minColIndex])

        if(matrix[minRowIndex][minColIndex] > 1):
            print("Devided")
            divider = matrix[minRowIndex][minColIndex]
            for i in range(0, len(matrix[0])):
                matrix[minRowIndex][i] /= divider


        for i in range(0, len(matrix)):
            if(i != minRowIndex):
                print("i: ", i)
                multiplier = matrix[i][minColIndex] / matrix[minRowIndex][minColIndex]
                for j in range(0, len(matrix[0])):
                    matrix[i][j] -= multiplier * matrix[minRowIndex][j]
    

        print(matrix)   

def main():
    #a,b,c = 0, 0, 4 # -20
    #       6, 0, 6 # -30
    a,b,c = 8, 10, 3
    #a, b, c = 9, 5, 6

    matrix =[[ 0,  2, -3,  0, -5,  0,  0,  0],
             [ a, -1,  1, -1, -1,  1,  0,  0], 
             [ b,  2,  4,  0,  0,  0,  1,  0],
             [ c,  0,  0,  1,  1,  0,  0,  1]]

    matrix2 =[[ 0,  3, -2, -6,  0,  0,  0,  0],
              [12,  4,  2,  0,  1,  1,  0,  0], 
              [ 8, -2,  2,  1, -1,  0,  1,  0],
              [ 5,  0,  0,  1,  2,  0,  0,  1]]
    
    #print(list(range(1, len(matrix))))
    simplex_method(matrix)


if __name__ == "__main__":
    main()
