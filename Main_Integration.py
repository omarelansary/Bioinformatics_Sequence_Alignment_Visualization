from GUIV1_4 import Ui_MainWindow
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow,QFileDialog,QMessageBox,QWidget
from PyQt5.QtGui import QPixmap
from mplwidget import MplWidget 
import qdarkstyle
import sys
from pathlib import Path
import os
from Global_Int_Version import *
from MSA import *
from Local_Int_Version import*
import subprocess
from Bio import SeqIO




class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        #Variales
        self.gap=0
        self.match=0
        self.mismatch=0
        self.All_Seq_Arr=[]

        # Show and Hide Initially
        self.ui.Browse_inputtype_pushButton.hide()
        self.ui.BrowsedFile_inputtype_label.hide()
        self.ui.groupBox_2.hide() # Enter text
        self.ui.GlobalAlignment_radioButton.hide()
        self.ui.LocalAlignment_radioButton.hide()
        self.ui.GapScore_label.hide()
        self.ui.GapScore_lineEdit.hide()
        self.ui.MatchScore_label.hide()
        self.ui.MatchScore_lineEdit.hide()
        self.ui.MismatchScore_label.hide()
        self.ui.MismatchScore_lineEdit.hide()
        
        self.ui.SumOFpairs_Score_label.hide()
        self.ui.PercentIdentify_Score_label.hide()
        self.ui.Mutualinformation_Score_label.hide()
        self.ui.NormalizedMutualInformation_Score_Label.hide()


        #Show and Hide when annotation
        self.ui.TotalScore_label.hide()#####
        self.ui.seq1_label.hide()
        self.ui.seq1_spinBox.hide()
        self.ui.seq2_label.hide()
        self.ui.seq2_spinBox.hide()


        self.ui.browse_inputtype_radioButton.clicked.connect(self.showBrowse)
        self.ui.textInput_inputtype_radioButton.clicked.connect(self.showEnterText)
        self.ui.pairWiseAlignment_radioButton.clicked.connect(lambda: self.showAlignmentType(0))
        self.ui.MultipleSeqAlignment_radioButton.clicked.connect(lambda: self.showAlignmentType(1))
        self.ui.ok_Cancel_buttonBox.accepted.connect(lambda: self.EnteredText(0))
        self.ui.ok_Cancel_buttonBox.rejected.connect(lambda: self.EnteredText(1))
        self.ui.GapScore_lineEdit.returnPressed.connect(lambda: self.ScoresLineEdit(0))
        self.ui.MatchScore_lineEdit.returnPressed.connect(lambda: self.ScoresLineEdit(1))
        self.ui.MismatchScore_lineEdit.returnPressed.connect(lambda: self.ScoresLineEdit(2))
        self.ui.Browse_inputtype_pushButton.clicked.connect(self.Browsefile)

       
        self.ui.actionShow_Input.triggered.connect(lambda: self.showandhideInput(0))
        self.ui.actionHide_Input.triggered.connect(lambda: self.showandhideInput(1))
        self.ui.actionRun_Alignment.triggered.connect(self.RunAlignmentCode)



    def showBrowse(self):
        self.ui.Browse_inputtype_pushButton.show()
        self.ui.BrowsedFile_inputtype_label.show()
        self.ui.groupBox_2.hide() # Enter text

    def showEnterText(self):
        self.ui.groupBox_2.show() # Enter text 
        self.ui.Browse_inputtype_pushButton.hide()
        self.ui.BrowsedFile_inputtype_label.hide()    

    def showAlignmentType(self,x):
        if (x==0):
            self.ui.TotalScore_label.show()#####
            self.ui.seq1_label.show()
            self.ui.seq1_spinBox.show()
            self.ui.seq2_label.show()
            self.ui.seq2_spinBox.show()

            self.ui.GlobalAlignment_radioButton.show()
            self.ui.LocalAlignment_radioButton.show()
            self.ui.GapScore_label.show()
            self.ui.GapScore_lineEdit.show()
            self.ui.MatchScore_label.show()
            self.ui.MatchScore_lineEdit.show()
            self.ui.MismatchScore_label.show()
            self.ui.MismatchScore_lineEdit.show()

            self.ui.SumOFpairs_Score_label.hide()
            self.ui.PercentIdentify_Score_label.hide()
            self.ui.Mutualinformation_Score_label.hide()
            self.ui.NormalizedMutualInformation_Score_Label.hide()
        else:
            self.ui.TotalScore_label.hide()#####
            self.ui.seq1_label.hide()
            self.ui.seq1_spinBox.hide()
            self.ui.seq2_label.hide()
            self.ui.seq2_spinBox.hide()

            self.ui.GlobalAlignment_radioButton.hide()
            self.ui.LocalAlignment_radioButton.hide()
            self.ui.GapScore_label.hide()
            self.ui.GapScore_lineEdit.hide()
            self.ui.MatchScore_label.hide()
            self.ui.MatchScore_lineEdit.hide()
            self.ui.MismatchScore_label.hide()
            self.ui.MismatchScore_lineEdit.hide()

            self.ui.SumOFpairs_Score_label.show()
            self.ui.PercentIdentify_Score_label.show()
            self.ui.Mutualinformation_Score_label.show()
            self.ui.NormalizedMutualInformation_Score_Label.show()

    def EnteredText(self,x):
        if (x==0):
            # Program to show various ways to read and
            # write data in a file.
            
            inputstring=str(self.ui.Entertext_textEdit.toPlainText()).upper()
            # sequences=[]
            # if (inputstring.find(">")!=-1):
            #     Firstline=int(inputstring.find(">"))
            #     newlineIndex=int(inputstring.find("\n"))
            #     Acsnnumber=inputstring[Firstline:(newlineIndex)]
            #     sequences=inputstring[(newlineIndex+1):].upper()

            var=inputstring.split(" ")
            file1 = open("myfile.txt","w")
            lines=[]
            self.All_Seq_Arr=[]#remove every thing on the array
            for i in range(len(var)):
                # print(str(i),var[i])
                # if (var[i].find(">")!=-1):
                #     print("found")
                #     lines.append(var[i])
                # else:
                #     print("not found")       
                lines.append(">"+str(i+1)+" Out\n")
                string=var[i].replace('\n', '')
                self.All_Seq_Arr.append(string)
                Numberoflines=int(len(string)//70) #max char in a line is 80, so we can make it 70 accoring to fasta file sent
                for j in range(Numberoflines):
                    lines.append(string[int(70*(j)):int(70*(j+1))]+"\n")
                lines.append(string[int(Numberoflines*70):int(len(string)+1)]+"\n")
                

            # \n is placed to indicate EOL (End of Line)
            file1.writelines(lines)
            file1.close() #to change file access modes

            # If file exists, delete it.
            if os.path.isfile("myfile.fasta"):
                os.remove("myfile.fasta")
            else:
                pass

            p = Path('myfile.txt')
            p.rename(p.with_suffix('.fasta'))

        else:
            self.ui.Entertext_textEdit.clear()

    def ScoresLineEdit(self,x):
        if (x==0):
            if (self.ui.GapScore_lineEdit.text().find("-")!=-1):
                val=self.ui.GapScore_lineEdit.text()
                len(val)
                self.gap=-int(val[1:])
                print(self.gap)
            else:
                self.gap=int(self.ui.GapScore_lineEdit.text())
                print(self.gap)
            
            
                
        elif(x==1):
            # self.match=int(str(self.ui.MatchScore_lineEdit.text()))
            if (self.ui.MatchScore_lineEdit.text().find("-")!=-1):
                val=self.ui.MatchScore_lineEdit.text()
                len(val)
                self.match=-int(val[1:])
                print(self.match)
            else:
                self.match=int(self.ui.MatchScore_lineEdit.text())
                print(self.match)
        else:
            # self.mismatch=int(str(self.ui.MismatchScore_lineEdit.text()))
            if (self.ui.MismatchScore_lineEdit.text().find("-")!=-1):
                val=self.ui.MismatchScore_lineEdit.text()
                len(val)
                self.mismatch=-int(val[1:])
                print(self.mismatch)
            else:
                self.mismatch=int(self.ui.MismatchScore_lineEdit.text())
                print(self.mismatch)
            

    def showandhideInput(self,x):
        if(x==0):
           self.ui.InputWindow_groupBox.show()
        else:
            self.ui.InputWindow_groupBox.hide()    

    def Browsefile(self):
        self.path = QFileDialog.getOpenFileName(self, 'Open a file', '') #open browse window
        if self.path != ('', ''):
            self.data = self.path[0]
            head, tail = os.path.split(self.data)
            if (tail.find(".fasta")==-1):
                self.data=[]
                self.messagebox("Please browse a .fasta file to Continue.")
            else:
                self.ui.BrowsedFile_inputtype_label.setText(tail)    
                records = list(SeqIO.parse(self.data, "fasta"))
                self.All_Seq_Arr=[]#clear and add
                for i in range(len(records)):
                    self.All_Seq_Arr.append(records[i].seq)

    def RunAlignmentCode(self):
        print("Running")
        
        if (len(self.All_Seq_Arr)!=0):
            seq1index=int(self.ui.seq1_spinBox.value())
            seq2index=int(self.ui.seq2_spinBox.value())
            if (seq1index>=len(self.All_Seq_Arr))or(seq2index>=len(self.All_Seq_Arr)):
                outputmsg="Error sequences index out of range ,the range is from 0 to: "+str(len(self.All_Seq_Arr)-1)
                self.messagebox(outputmsg)
            else:
                if (self.ui.GlobalAlignment_radioButton.isChecked()):
                    Glob_path,Total_Scores=Global_Alignment(self.All_Seq_Arr[seq1index],self.All_Seq_Arr[seq2index],self.match,self.mismatch,self.gap)
                    self.ui.AlignmentDisplay_label.setPixmap(QPixmap("Global_Align.png"))   
                    Output_TotalScore= 'Total Score: '+str(Total_Scores)   
                    self.ui.TotalScore_label.setText(Output_TotalScore)
                    print("Finished Global Align")
                    
                if (self.ui.LocalAlignment_radioButton.isChecked()):
                    Local_path,Total_Scores=Local_Alignment(self.All_Seq_Arr[seq1index],self.All_Seq_Arr[seq2index],self.match,self.mismatch,self.gap)
                    self.ui.AlignmentDisplay_label.setPixmap(QPixmap("Local_Align.png")) 
                    Output_TotalScore= 'Total Score: '+str(Total_Scores)   
                    self.ui.TotalScore_label.setText(Output_TotalScore)          
                    print("Finished Local Align")
            if (self.ui.MultipleSeqAlignment_radioButton.isChecked()):    
                # MSA
                if (self.ui.textInput_inputtype_radioButton.isChecked()):
                    MSA_path,percent_Identity,sum_of_pairs,MutualInfo,NormalizedMutualInfo=MSA("myfile.fasta")
                    self.ui.AlignmentDisplay_label.setPixmap(QPixmap("MultipleSeq_Align.png"))
                    print("Finished text input MSA")    
                elif(self.ui.browse_inputtype_radioButton.isChecked()):
                    MSA_path,percent_Identity,sum_of_pairs,MutualInfo,NormalizedMutualInfo=MSA(self.data)
                    self.ui.AlignmentDisplay_label.setPixmap(QPixmap("MultipleSeq_Align.png"))
                    print("Finished browse input MSA")

                output_MutualInfo="Mutual Infromation: "+str(round(MutualInfo,2))
                output_NormalizedMutualInfo="Normalized Mutual Infromation: "+str(round(NormalizedMutualInfo,2))    
                output_Percent_Identity="Percent Identity: "+str(round(percent_Identity,2))
                output_Sum_Of_pairs="Sum Of Pairs: "+str(sum_of_pairs)


                self.ui.PercentIdentify_Score_label.setText(output_Percent_Identity)
                self.ui.SumOFpairs_Score_label.setText(output_Sum_Of_pairs)
                self.ui.Mutualinformation_Score_label.setText(output_MutualInfo)
                self.ui.NormalizedMutualInformation_Score_Label.setText(output_NormalizedMutualInfo)    
                
            self.ui.BetterView_Alignment_label.setOpenExternalLinks(True)
            linkTemplate='<a href={0}>{1}</a>' 
            #put plot(fig) in each file Global,Local,MSA to return the output in temp-plot    
            #Please change this link to your own local host
            self.ui.BetterView_Alignment_label.setText(linkTemplate.format('file:///C:/Users/user/Desktop/Final_Integration/temp-plot.html','Better View click here'))
        else:
            self.messagebox("Enter or browse Sequences First.")

    def messagebox(self,x):    
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText(x)
        msg.setWindowTitle('ERROR')
        msg.exec_()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyleSheet(qdarkstyle.load_stylesheet())
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())      

