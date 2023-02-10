import numpy as np
import numpy
import plotly
import plotly.graph_objects as go
from plotly.offline import *
import plotly.figure_factory as ff
# match=1
# mismatch=-1
# gap=-2
dic={'A': 230, 'B': 20, 'C': 300, 'D': 40, 'E': 50, 'F': 60, 'G': 270, 'H': 80, 'I': 90, 'J':100, 'K': 110, 'L': 120, 'M': 130, 'N': 140, 'O': 150, 'P': 160, 'Q': 170, 'R': 180, 'S': 190, 'T': 150, 'U': 210, 'V': 220, 'W': 230, 'X': 235, 'Y': 240, 'Z': 245,'_':200,'Total':-1,'Scores':-1}
# dic={'A': 200, 'B': 20, 'C': 30, 'D': 40, 'E': 50, 'F': 60, 'G': 90, 'H': 80, 'I': 90, 'J':100, 'K': 110, 'L': 120, 'M': 130, 'N': 140, 'O': 150, 'P': 160, 'Q': 170, 'R': 180, 'S': 190, 'T': 200, 'U': 210, 'V': 220, 'W': 230, 'X': 235, 'Y': 240, 'Z': 245,'_':100,'Total':-1,'Scores':-1}
print(len(dic)) 
def checker_index_negative(i,j):
        if i<0 or j<0:
            return -1
def Compute_Scores(a,b,match,mismatch,gap):
    scores_list=[]
    colors_seq1=[]
    colors_seq2=[]
    for x, y in zip(a, b):
        col1=int(dic[x])
        col2=int(dic[y])
        colors_seq1.append(col1)
        colors_seq2.append(col2)
        if x == y:
            scores_list.append(match)
        elif x!=y:
            if x=='_' or y=='_':
                scores_list.append(gap)
            else:
                scores_list.append(mismatch)
    tot_score=np.sum(scores_list)
    scores_list.append(tot_score)
    # colors_seq1.append(dic['Total'])
    # colors_seq2.append(dic['Scores'])
    Matrix_colors=[
    
    colors_seq1,
    colors_seq2,
    ]
    # Matrix_colors=np.row_stack((colors_seq2,colors_seq1))
    
    # print("after matrix cols:")
    # print(Matrix_colors)
    return scores_list,Matrix_colors
def Global_Alignment(seq1,seq2,match,mismatch,gap):

    matrix=numpy.zeros((len(seq1)+1,len(seq2)+1))

    #initializations
    arr=numpy.zeros((3))

    for i in range(len(seq2)+1):
        matrix[0][i]=i*gap

    for i in range(len(seq1)+1):
        matrix[i][0]=i*gap

    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if (seq1[i-1]==seq2[j-1]):
                arr[0]=matrix[i-1][j-1]+match
            else:
                arr[0]=matrix[i-1][j-1]+mismatch
            arr[1]=matrix[i-1][j]+gap
            arr[2]=matrix[i][j-1]+gap
            matrix[i][j]=numpy.max(arr)        
    # print(matrix)

    seq2=[i for i in seq2]
    seq1=[i for i in seq1]

# This function gets the shortest path from the start to the goal
# def shortest_path(matrix):
    # Create list of neighbours
    list_neighbours = []
    # Create list of tuples for the indices of elements of path taken(included start and goal indices)
    list_tuples = []
    # [(i,j),(i2,j2)]
    # Acts as a flag when goal is found 
    goalfound = 0
    # Get the start indices of your matrix
    max = -99999
    i_start = len(seq1)
    j_start = len(seq2)
    # print("My start row: ", i_start)
    # print("My start col: ", j_start)
    # for i in range(matrix.shape[0]):
    #     for j in range(matrix.shape[1]):
    #         if (matrix[i][j]>max):
    #             max=matrix[i][j]
    #             #Assign start and end indices to i&j    
    #             i_start=i
    #             j_start=j
    # arr.append([i,j]) 
    # print("goal indices:",i,j)
    i = i_start
    j = j_start
    upperleftstart = matrix[i - 1][j - 1]
    leftstart = matrix[i][j - 1]
    upperstart = matrix[i - 1][j]
    if seq1[i-1]==seq2[j-1]:
        upperleftstart=upperleftstart+match
    elif seq1[i-1]!=seq2[j-1]:
        upperleftstart=upperleftstart+mismatch
    upperstart=upperstart+gap
    leftstart=leftstart+gap
    # Comparison to get maximum
    if upperleftstart >= leftstart:
        max1ini = upperleftstart
    else:
        max1ini = leftstart
    if max1ini >= upperstart:
        max2ini = max1ini
    else:
        max2ini = upperstart
    # max1ini=np.max(upperleftstart,upperstart)
    # max2ini=np.max(max1ini,leftstart)
    if max2ini == upperleftstart:
        # Make a tuple to add to list
        mytuple = (i, j, 'upperleft')
        # Append it to the list
        list_tuples.append(mytuple)
    elif max2ini == leftstart:
        # Make a tuple to add to list
        mytuple = (i, j, 'left')
        # Append it to the list
        list_tuples.append(mytuple)
    elif max2ini == upperstart:
        # Make a tuple to add to list
        mytuple = (i, j, 'upper')
        # Append it to the list
        list_tuples.append(mytuple)
    # print('i=i_start= ', i)
    # print('j=j_start= ', j)

    # Loop while goal not found
    while (i != 0 and j != 0):
        # If goal found break the while loop
        # if(matrix[i][j]==0):
        #     goalfound=1
        #     break
        list_neighbours = []
        # Get neighbours and check for each if the index is negative (out of boundaries)
        checker = checker_index_negative(i - 1, j)  # Upper Gap
        if (checker != -1):
            upper = matrix[i - 1][j]
            upper=upper+gap
            list_neighbours.append(upper)
        else:
            upper = -999999999999

        #
        checker = checker_index_negative(i, j - 1)  # Left gap
        if (checker != -1):
            left = matrix[i][j - 1]
            left=left+gap
            list_neighbours.append(left)
        else:
            left = -999999999999

        checker = checker_index_negative(i - 1, j - 1)  # Upper left (match or mismatch)
        if (checker != -1):
            upperleft = matrix[i - 1][j - 1]
            #check if match or mismatch
            if seq1[i-1]==seq2[j-1]:
                upperleft=upperleft+match
            elif seq1[i-1]!=seq2[j-1]:
                upperleft=upperleft+mismatch
            list_neighbours.append(upperleft)
        else:
            upperleft = -999999999999

        # Loop over list of neighbours and execludes the ones(obstacles)
        # for k in range(len(list_neighbours)):
        #  if list_neighbours[k]!=1:
        #     list_neighbours_2.append(list_neighbours[k])           
        max_val = -999999999999
        # Get the minimium the neighbours to all neighbours except the ones(obstacles)
        max_val = np.max(list_neighbours)
        # print("Mymax=", max_val)
        
        if (upperleft == max_val):  # Upper left 
            i = i - 1
            j = j - 1
            # =============================================
            upperleftstart = matrix[i - 1][j - 1]
            leftstart = matrix[i][j - 1]
            upperstart = matrix[i - 1][j]
            if seq1[i-1]==seq2[j-1]:
               upperleftstart=upperleftstart+match
            elif seq1[i-1]!=seq2[j-1]:
                upperleftstart=upperleftstart+mismatch
            upperstart=upperstart+gap
            leftstart=leftstart+gap
            if upperleftstart >= leftstart:
                max1ini = upperleftstart
            else:
                max1ini = leftstart
            if max1ini >= upperstart:
                max2ini = max1ini
            else:
                max2ini = upperstart
            # max1ini=np.max(upperleftstart,upperstart)
            # max2ini=np.max(max1ini,leftstart)
            if max2ini == upperleftstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'upperleft')
                # Append it to the list
                list_tuples.append(mytuple)
            elif max2ini == leftstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'left')
                # Append it to the list
                list_tuples.append(mytuple)
            elif max2ini == upperstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'upper')
                # Append it to the list
                list_tuples.append(mytuple)
            ###33333333333333333333333333333333333333333

        # The upper is the max 
        # if(upper==max_val or upper==0):

        # The upperleft is the max 
        elif (upper == max_val):  # Upper left 
            i = i - 1
            j = j

            upperleftstart = matrix[i - 1][j - 1]
            leftstart = matrix[i][j - 1]
            upperstart = matrix[i - 1][j]
            if seq1[i-1]==seq2[j-1]:
               upperleftstart=upperleftstart+match
            elif seq1[i-1]!=seq2[j-1]:
                upperleftstart=upperleftstart+mismatch
            upperstart=upperstart+gap
            leftstart=leftstart+gap
            if upperleftstart >= leftstart:
                max1ini = upperleftstart
            else:
                max1ini = leftstart
            if max1ini >= upperstart:
                max2ini = max1ini
            else:
                max2ini = upperstart
            # max1ini=np.max(upperleftstart,upperstart)
            # max2ini=np.max(max1ini,leftstart)
            if max2ini == upperleftstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'upperleft')
                # Append it to the list
                list_tuples.append(mytuple)
            elif max2ini == leftstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'left')
                # Append it to the list
                list_tuples.append(mytuple)
            elif max2ini == upperstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'upper')
                # Append it to the list
                list_tuples.append(mytuple)
            ###33333333333333333333333333333333333333333


        

        elif (left == max_val):  # left 
            i = i
            j = j - 1
            
            upperleftstart = matrix[i - 1][j - 1]
            leftstart = matrix[i][j - 1]
            upperstart = matrix[i - 1][j]
            if seq1[i-1]==seq2[j-1]:
               upperleftstart=upperleftstart+match
            elif seq1[i-1]!=seq2[j-1]:
                upperleftstart=upperleftstart+mismatch
            upperstart=upperstart+gap
            leftstart=leftstart+gap
            if upperleftstart >= leftstart:
                max1ini = upperleftstart
            else:
                max1ini = leftstart
            if max1ini >= upperstart:
                max2ini = max1ini
            else:
                max2ini = upperstart
            # max1ini=np.max(upperleftstart,upperstart)
            # max2ini=np.max(max1ini,leftstart)
            if max2ini == upperleftstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'upperleft')
                # Append it to the list
                list_tuples.append(mytuple)
            elif max2ini == leftstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'left')
                # Append it to the list
                list_tuples.append(mytuple)
            elif max2ini == upperstart:
                # Make a tuple to add to list
                mytuple = (i, j, 'upper')
                # Append it to the list
                list_tuples.append(mytuple)
                ###33333333333333333333333333333333333333333

    seq1_list = []
    seq2_list = []
    val1 = 'u'
    val2 = 'u'
    # print(len(list_tuples))
    # print(list_tuples[0][0] - 1)
    # print(list_tuples)
    for i in range(len(list_tuples)):
        if list_tuples[i][2] == 'upperleft':  # list_tuples[i][0]==list_tuples[i][1]:
            val1 = seq1[(list_tuples[i][0]) - 1]
            val2 = seq2[(list_tuples[i][1]) - 1]
            # print("index1: ", (list_tuples[i][0]) - 1)
            # print("index2: ", (list_tuples[i][1]) - 1)
            seq1_list.append(val1)
            seq2_list.append(val2)
        elif list_tuples[i][
            2] == 'left':  # ((list_tuples[i][0])>(list_tuples[i][1])):#seq1 is greater index so gap in seq2
            val1 = '_'
            val2 = seq2[(list_tuples[i][1]) - 1]
            seq1_list.append(val1)
            seq2_list.append(val2)
            # print("index1: ",(list_tuples[i][0])-1)
            # print("index2: ", (list_tuples[i][1]) - 1)
        elif list_tuples[i][2] == 'upper':  # ((list_tuples[i][0])<(list_tuples[i][1])):
            val1 = seq1[(list_tuples[i][0]) - 1]
            val2 = '_'
            seq1_list.append(val1)
            seq2_list.append(val2)
            # print("index1: ", (list_tuples[i][0]) - 1)
            # print("index2: ",(list_tuples[i][1])-1)
    seq1f = []
    seq2f = []
    seq2f = seq2_list[:len(seq2_list) - 1]
    seq1f = seq1_list[:len(seq1_list) - 1]
    # print('seq2_list: ', seq2f[::-1])
    # print('seq1_list: ', seq1f[::-1])
    seq2f=seq2f[::-1]
    seq1f=seq1f[::-1]
    seq_2f_string = "".join(seq2f)
    seq_1f_string= "".join(seq1f)
    print('seq2_list: ',seq_2f_string)
    print('seq1_list: ',seq_1f_string)
    # print("Seq2f:")
    # print(seq2f)
    # print("Seq1f:")
    # print(seq1f)
    #===============================Original
    ScoresList,Matrix_colors=Compute_Scores(seq1f,seq2f,match,mismatch,gap)
    Tot_scores_val=ScoresList[-1]
    # ''.join(ScoresList)
    ScoresList=[str(i) for i in ScoresList]
    print("Matrix after call:")
    print(Matrix_colors)
    # seq2f.append("Total")
    # seq1f.append("Scores")
    fig = go.Figure(data=go.Heatmap(
                    z=Matrix_colors
                                   ,
                    text=[
                          seq1f, 
                          seq2f,                                                                           
                          ],
                    xgap = 5,
                    ygap = 5,
                    texttemplate="%{text}",
                    textfont={"size":3}
                    , coloraxis=None, colorbar=None)
                    ,layout=go.Layout(height=700, width=700))

    # z=Matrix_colors

    # text=[
    #                       seq1f, 
    #                       seq2f,                                                                           
    #                       ScoresList
    #                       ]
    # fig = ff.create_annotated_heatmap(z, annotation_text=text, colorscale='Viridis',
    #                               hoverinfo='z')
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_coloraxes(showscale=False)
    

    mattrx=[
          seq2f, 
          seq1f,                                 
          ScoresList
    ]
    print("Matrix of Seqs:")
    print(mattrx)
    # fig.write_image("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # iplot(fig)
    fig.write_image("Global_Align.png")
    # plotly.offline.plot(fig)
    # fig.save("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    #fig.show()
    
    #==============================================
    # Return list of tuples for indices of my path
    list_tuples=list_tuples[::-1]
    print(list_tuples[:len(list_tuples) - 1])
    path=plot(fig)
    #=========End original
    #===================================================
    # #=================Part One
    # print("After Align Length:",len(seq1f))

    # ScoresList,Matrix_colors=Compute_Scores(seq1f,seq2f,match,mismatch,gap)
    # # ''.join(ScoresList)
    # ScoresList=[str(i) for i in ScoresList]
    # print("Matrix after call:")
    # print(Matrix_colors)
    # # seq2f.append("Total")
    # # seq1f.append("Scores")

    # seq1f1=seq1f[0:200]
    # print("seq[0:200]")
    # print(seq1f1)
    # seq2f1=seq2f[0:200]
    # ScoresList1=ScoresList[0:200]
    # Matrix1=[
    #     [1,2,3],
    #     [3,4,5],
    #     [6,7,8]
    # ]
    # Matrix1[0]=Matrix_colors[0][0:200]
    # Matrix1[1]=Matrix_colors[1][0:200]
    # Matrix1[2]=Matrix_colors[2][0:200]



    # fig = go.Figure(data=go.Heatmap(
    #                 z=Matrix1
    #                                ,
    #                 text=[
    #                       seq1f1, 
    #                       seq2f1,                                                                           
    #                       ScoresList1
    #                       ],
    #                 xgap = 5,
    #                 ygap = 5,
    #                 texttemplate="%{text}",
    #                 textfont={"size":5}
    #                 , coloraxis=None, colorbar=None)
    #                 ,layout=go.Layout(height=700, width=700))

    # fig.update_xaxes(visible=False)
    # fig.update_yaxes(visible=False)
    # fig.update_coloraxes(showscale=False)
    

    # # mattrx=[
    # #       seq2f, 
    # #       seq1f,                                 
    # #       ScoresList
    # # ]
    # # print("Matrix of Seqs:")
    # # print(mattrx)
    # # fig.write_image("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # # iplot(fig)
    # fig.write_image("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align1.png")
    # # plotly.offline.plot(fig)
    # # fig.save("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # # fig.show()
    # #==============================================
    # # Return list of tuples for indices of my path
    # list_tuples=list_tuples[::-1]
    # print(list_tuples[:len(list_tuples) - 1])
    # path1="E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align1.png"
    # #========================End Part one
    # #==================Part two
    # seq1f2=seq1f[200:400]
    # print("seq2[200:400]")
    # print(seq1f2)
    # seq2f2=seq2f[200:400]
    # ScoresList2=ScoresList[200:400]
    # Matrix2=[
    #     [1,2,3],
    #     [3,4,5],
    #     [6,7,8]
    # ]
    # Matrix2[0]=Matrix_colors[0][200:400]
    # Matrix2[1]=Matrix_colors[1][200:400]
    # Matrix2[2]=Matrix_colors[2][200:400]



    # fig = go.Figure(data=go.Heatmap(
    #                 z=Matrix2
    #                                ,
    #                 text=[
    #                       seq1f2, 
    #                       seq2f2,                                                                           
    #                       ScoresList2
    #                       ],
    #                 xgap = 5,
    #                 ygap = 5,
    #                 texttemplate="%{text}",
    #                 textfont={"size":5}
    #                 , coloraxis=None, colorbar=None)
    #                 ,layout=go.Layout(height=700, width=700))

    # fig.update_xaxes(visible=False)
    # fig.update_yaxes(visible=False)
    # fig.update_coloraxes(showscale=False)
    

    # # mattrx=[
    # #       seq2f, 
    # #       seq1f,                                 
    # #       ScoresList
    # # ]
    # # print("Matrix of Seqs:")
    # # print(mattrx)
    # # fig.write_image("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # # iplot(fig)
    # fig.write_image("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align2.png")
    # # plotly.offline.plot(fig)
    # # fig.save("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # # fig.show()
    # #==============================================
    # # Return list of tuples for indices of my path
    # list_tuples=list_tuples[::-1]
    # print(list_tuples[:len(list_tuples) - 1])
    # path2="E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align2.png"
    # #end part 2 

    print(path)
    return path,Tot_scores_val


#===================================================================================================================================

# seq1="ACGT"
# seq2="A_GC"
# seq2=[i for i in seq2]
# seq1=[i for i in seq1]

#=====================================Back up for fig
# ScoresList,Matrix_colors=Compute_Scores(seq1f,seq2f,match,mismatch,gap)
    # # ''.join(ScoresList)
    # ScoresList=[str(i) for i in ScoresList]
    # print("Matrix after call:")
    # print(Matrix_colors)
    # seq2f.append("Total")
    # seq1f.append("Scores")
    # fig = go.Figure(data=go.Heatmap(
    #                 z=Matrix_colors
    #                                ,
    #                 text=[
    #                       seq1f, 
    #                       seq2f,                                                                           
    #                       ScoresList
    #                       ],
    #                 xgap = 5,
    #                 ygap = 5,
    #                 texttemplate="%{text}",
    #                 textfont={"size":5}
    #                 , coloraxis=None, colorbar=None)
    #                 ,layout=go.Layout(height=700, width=700))

    # # z=Matrix_colors

    # # text=[
    # #                       seq1f, 
    # #                       seq2f,                                                                           
    # #                       ScoresList
    # #                       ]
    # # fig = ff.create_annotated_heatmap(z, annotation_text=text, colorscale='Viridis',
    # #                               hoverinfo='z')
    # fig.update_xaxes(visible=False)
    # fig.update_yaxes(visible=False)
    # fig.update_coloraxes(showscale=False)
    

    # mattrx=[
    #       seq2f, 
    #       seq1f,                                 
    #       ScoresList
    # ]
    # print("Matrix of Seqs:")
    # print(mattrx)
    # # fig.write_image("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # # iplot(fig)
    # fig.write_image("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # # plotly.offline.plot(fig)
    # # fig.save("E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png")
    # fig.show()
    # #==============================================
    # # Return list of tuples for indices of my path
    # list_tuples=list_tuples[::-1]
    # print(list_tuples[:len(list_tuples) - 1])
    # path="E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Final_Integration\\Global_Align.png"

#=======================================================


# def Alignment_Visualization(seq1,seq2,scores_list):
    
#     # Visualization_Matrix=np.array([[seq1],[seq2],[scores_list]],dtype=object)
#     Visualization_Matrix=[[seq1],[seq2],[scores_list]]
#     return Visualization_Matrix

# print("Scores:")
# print(ScoresList)

# matrix=Alignment_Visualization(seq1,seq2,ScoresList)
# print("Matrix:")
# print(matrix)
# print("Seq1colors:")
# print(colors_seq1)
# print("Seq2colors")
# print(colors_seq2)
# print("list zefftt")
# print(Matrix_colors)
# print(type(Matrix_colors))
#================================



# import plotly.express as px

# # z = [[.1, .3, .5, .7, .9],
# #      [1, .8, .6, .4, .2],
# #      [.2, 0, .5, .7, .9],
# #      [.9, .8, .4, .2, 0],
# #      [.3, .4, .5, .7, 1]]
# # z=[
# # [1,3,4,5,7,9],

# # ]

# fig = px.imshow(Matrix_colors, text_auto=True)
# fig.show()