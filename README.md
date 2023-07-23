# Dislocation Density VUMAT Subroutine for ABAQUS

This is a simple VUMAT subroutine written in FORTRAN to implement a Dislocation Density Constitutive Model in ABAQUS.
The subroutine is quite simple that only a few basic FORTRAN knowledge needed for the VUMAT, and some rules about program which ABAQUS made is also essential for understanding. These program rules, which mostly about the input and output parameters of subroutine or the information of how the parameters are set, could be find in the ABAQUS Documentation.

**To say the least, if you didn't want to read the code and understand how it works, The need of these knowledge could be ignored.** You just need to load the *.for file to the ABAQUS with a VUMAT environment.

It is needed to mentioned that this VUMAT subroutine isn't perfect now. The convergence of some governing equation used in Dislocation Density Constitutive Model can not be guaranteed. Some details about the constitutive model could be find in literature with DOI: https://doi.org/10.1016/j.ijmecsci.2018.08.005 (I'm not the author of this literature. And appreciate to its authors, my code was written according to their works).

I just put my work here. And I think it is a good start could be utilized if you were new who doing some research about the Dislocation Density Constitutive Model too. I suffered from few study material in VUMAT subroutine. Hope this could be helpful and welcome to communicate with my email-(a87669289@163.com) if there is any problem.
