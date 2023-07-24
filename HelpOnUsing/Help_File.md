# The example of Using dislocation.for

#### It need to be mentioned that my ABAQUS version is 6.14-4, although the version might doesn't matter. 

### 1. How to load the dislocation.for?

- #### The Fortran Compiler

Firstly, you should make sure that your ABAQUS already had a "VUMAT environment". You can run a "system configuration checks" in the **windows command line** by the code `abaqus info=system` to check if your ABAQUS has a Fortran compiler which is needed for our subroutine. 

```
# ABAQUS without Fortran Compiler
>>> abaqus info=system
...
Fortran Compiler:     Unable to find a Fortran compiler on this system.  If
                      Intel Fortran is installed on this system, please load
                      ifortvars.bat before running Abaqus
...(some othe information)

# ABAQUS with Fortran Compiler
>>> abaqus info=system
...
Fortran Compiler:     Intel Fortran Compiler 14.0
...(some othe information)
```

You can find how to install a Fortran Compiler and how to link it with the ABAQUS in some other tutorial. It is a quite complex thing, so we assume that your ABAQUS already had a Fortran Compiler and do not introduce how to install it here.



- #### Load the dislocation.for

The following picture shows the detail operation.

<img src="\imgs\01.png" alt="image01" style="zoom: 80%;" />





### 2. Define the Parameters

- #### The Parameters' meaning and the Value of them

According to another literature (http://dx.doi.org/10.1016/j.actamat.2015.06.054), the parameters of the constitutive model are the following picture.

<img src="\imgs\02.png" alt="image02" style="zoom: 80%;" />

These parameters values are the values which there works use. You can adjust them according to your task.

And the Tuned parameters is not related to the property of the material, while the Other parameters is related to the property of  the material. You can just use the values above to do a test.

My code uses 22 parameters, you can find the meaning of each of them in the source code. But there are only 20 parameters in the picture, that because I add the E_target (Young Modulus), the Mu_target (Poisson's Ratio) and the S1_target (Yield Stress), while I remove the G (Shear Modulus) in the picture parameters (because it can be calculated with the Young Modulus and the Poisson's Ratio). So in My code there are 22 parameters.

To make it convenience, I transfer these parameters in the picture to the text so that you can copy the value of them. 

```fortran
E_target = 200000.0
Mu_target = 0.22
S1_target = 792.0
Alpha_star = 0.154
Beta_star = 0.078
K_c = 18.6
K_w = 32.8
N_c = 89.8
N_w = 90.3
M_star = 60.8
Burger = 2.48e-7
Taylor = 3.06
Eta = 0.25
VolFrac_zero = 0.25
VolFrac_infin = 0.06
Gamma0 = 1.0e7
K_zero = 100.0
K_infin = 1.0
Peeq_wave = 3.2
Beta = 0.26
Rho_c0 = 2.5e7
Rho_w0 = 5.0e7
```



- #### Define the parameters in the ABAQUS

To define the parameters which we use in out subroutine, you show define a material in the *Property* view firstly. It's as the same as how you define a material before. Our VUMAT subroutine is for ABAQUS/Explicit, so that define a ***Density*** material behavior is essential. Except the ***Density***, we have the other two material behavior which are the ***Depvar*** and the ***User Material***. They are all in the pull-down option under the *General* option. The detail operation the value of them is shown in the following picture.

<center>
    <img src="\imgs\03.png" alt="image03" style="zoom:50%;" />
    <img src="\imgs\04.png" alt="image04" style="zoom:50%;" />
</center>
<img src="\imgs\05.png" alt="image-05" style="zoom:67%;" />

***Depvar*** is the solution dependent variables we define in our VUMAT subroutine. And there are 7 variables I defined, so that the value of the ***Depvar*** is 7. And the meaning of each variable is appeared in the commentary within the source code. I also put them in the following table. These information can help you to check each of these variable while knowing their meaning during the post-processing.

| Variable |               Meaning               |
| :------: | :---------------------------------: |
|   SDV1   |            Yield Stress             |
|   SDV2   |                PEEQ                 |
|   SDV3   |      Dislocation Wall Density       |
|   SDV4   |      Dislocation Cell Density       |
|   SDV5   |      Dislocation Cell Diameter      |
|   SDV6   | Volume Fraction of Dislocation Wall |
|   SDV7   |      Total Dislocation Density      |

By the way, the order of the parameters in the ***User Material*** is important, make sure they are not in wrong order. And you may notice that some of these value are different with which in the literature. For example, the No.11 parameter is 2.48E-007 but it is 2.48E-010 in the literature. That's because we use centimeter as unit in our ABAQUS while the unit is meter in the literature, so that these are some conversion in there.



### 3. Test Model

When you are finishing to define the parameters, then you can use this Dislocation Density Constitutive Model actually. So I will quickly to  introduce a Test Model which I used to test our VUMAT subroutine and some of its effect.

This part is mainly to show the effect of our VUMAT subroutine. So to make the Help file briefly, some details of modeling will be ignored. And this is the last part of this Help file, don't be afraid, the test model is as simple as you can imagine.

- #### the geometry and the load

Our test model is a simple cantilever beam, its right side were fixed. Then we define a reference point **RP-1** and coupling it to the left side of the beam so that we can control the left side by the **RP-1**.

The load is applied on the **RP-1**. Actually, its a Boundary Condition in ABAQUS. we fixed the **RP-1** on all the dimension unless the UR3 which mean the rotate on the axis Z. The value of UR3 is modified with the value **6.28** in the analyse step which means we let the left side rotate 2 circle. So the cantilever beam will undergo torsion deformation which will create lots of shear stress, so the dislocation density will vary greatly.

<img src="\imgs\06.png" alt="image06" style="zoom:80%;" />



- #### Adding the SDV into the Field Output

**It need to be mentioned that adding the SDV into the Field Output is a important thing if you want to check the 7 solution dependent variables we mentioned before.** The detail of operation is shown in the following picture.

<img src="\imgs\07.png" alt="image07" style="zoom:80%;" />



- #### The effect

Mises:

<img src="\imgs\08.png" alt="image08" style="zoom:80%;" />



SDV5 (Cell Size):

<img src="\imgs\09.png" alt="image09" style="zoom:80%;" />



SDV7 (Dislocation Density):

<img src="\imgs\10.png" alt="image10" style="zoom:80%;" />



**The test_dislocation.cae and the job-1.odb has been uploaded on the same folder of this Help file. You can download it for learning.**

#### Thanks for Reading

