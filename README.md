Introduction
------------
![](https://github.com/SunXQlab/CCCbank/blob/main/CCCbank.png)

To facilitate practical applications of CCC inference tools, we
developed an integrated computational pipeline, CCCbank, that enables
flexible combinations between 16 LR inference methods and 13 LR
databases. Via CCCbank, each method can be equipped with one or more LR
databases or any other user-defined databases. CCCbank uses Seurat
object of scRNA-seq data as input for CCC inference. The output of
CCCbank includes inferred LR pairs and their scores, which can be used
for further analysis or visualization.

Installation
------------

CCCbank can be installed from the github.

    if(!require(devtools)){
      install.packages('devtools')
    }
    devtools::install_github('SunXQlab/CCCbank')

Usage
-----

Run CCC inference methods with their default LR prior database, for
example, CellChat.

    # Method 1:
    result <- CCCbank(ser = ser, method = 'CellChat', species = 'human')

    # Method 2:
    result <- RunCellChat(ser = ser, species = 'human')

If you want to run other methods like method 2 above, please refer to
the corresponding RunFunction in the table in ***Methods***.

Run CCC inference methods with one or more built-in LR prior databases.
Note that if using built-in LR prior databases, species is only ‘human’.

    # run CellChat method with the combination of CellPhoneDBLR and FANTOM5 
    # as LR prior database to infer CCC.
    result <- CCCbank(ser = ser, method = 'CellChat', species = 'human', 
                      databases = c('CellPhoneDBLR', 'FAMTOM5'), extension = FALSE)

    # run CellChat method with the combination of CellPhoneDBLR, FANTOM5 
    # and defualt LR prior database of CellChat as LR prior database
    result <- CCCbank(ser = ser, method = 'CellChat', species = 'human', 
                      databases = c('CellPhoneDBLR', 'FAMTOM5'), extension = TRUE)

Run CCC inference methods with user-defined databases.

    database <- data.frame(ligand = c('L1', 'L2', 'L3', 'L4'),
                           receptor = c('R1', 'R2', 'R1', 'R3'))
    # Method 1
    result <- CCCbank(ser = ser, method = 'CellChat', species = 'human', 
                      databases = database, extension = FALSE)

    # Method 2
    database <- data.frame(ligand = c('L1', 'L2', 'L3', 'L4'),
                           receptor = c('R1', 'R2', 'R1', 'R3'),
                           lr = c('L1_R1', 'L2_R2', 'L3_R1', 'L4_R3'))
    db_new <- ChangeCellChatDB(priorDatabase = database, 
                               extension = FALSE, species = 'human')
    result <- RunCellChat(ser = ser, species = 'human', priorDatabase = db_new)

If you want to run other methods with user-defined databases like method
2 above, please refer to the corresponding RunFunction and
ChangeDBFunction in the table in ***Methods***. And check the usage of
methods using ‘?RunFunction’, for example, ?RunCellChat. Note that if
the user-defined databases contain multi-subunit of ligands and
receptors, use ‘&’ to connect them. For example, if the ligand is A&B
and the receptor is D&W, lr is A&B\_D&W.

Methods
-------

There are **16 CCC inference methods** implemented in CCCbank R package:

<body>
<table border="1" cellspacing="1" width="800">
<thead>
<tr>
<th>
Methods
</th>
<th>
Version
</th>
<th>
Supported Species
</th>
<th>
RunFunction
</th>
<th>
ChangeDBFunction
</th>
</tr>
</thead>
<tbody>
<tr class="even-row">
<td>
<a href="https://github.com/ventolab/CellphoneDB">CellPhoneDB 2/3</a>
</td>
<td>
3.0.0
</td>
<td>
human
</td>
<td>
RunCellPhoneDB
</td>
<td>
ChangeCellPhoneDB
</td>
</tr>
<tr>
<td>
<a href="https://github.com/arc85/celltalker">CellTalker</a>
</td>
<td>
0.0.4.9000
</td>
<td>
human
</td>
<td>
RunCellTalker
</td>
<td>
ChangeCellTalkerDB
</td>
</tr>
<tr class="even-row">
<td>
<a href="https://github.com/msraredon/Connectome">Connectome</a>
</td>
<td>
1.0.1
</td>
<td>
human<br>mouse</br>rat</br>pig</br>
</td>
<td>
RunConnectome
</td>
<td>
ChangeConnectomeDB
</td>
</tr>
<tr>
<td>
<a href="https://github.com/soumelis-lab/ICELLNET">ICELLNET</a>
</td>
<td>
1.0.1
</td>
<td>
human
</td>
<td>
RunICELLNET
</td>
<td>
ChangeICELLNETDB
</td>
</tr>
<tr class="even-row">
<td>
<a href="https://github.com/asrhou/NATMI">NATMI</a>
</td>
<td>
——
</td>
<td>
human<br>mouse</br>(21 species in total)</br>
</td>
<td>
RunNATMI
</td>
<td>
ChangeNATMIDB
</td>
</tr>
<tr>
<td>
<a href="https://github.com/Coolgenome/iTALK">iTALK</a>
</td>
<td>
0.1.0
</td>
<td>
human
</td>
<td>
RuniTALK
</td>
<td>
ChangeiTALKDB
</td>
</tr>
<tr class="even-row">
<td>
<a href="https://github.com/JonETJakobsson/scConnect">scConnect</a>
</td>
<td>
1.0.3
</td>
<td>
human<br>mouse</br>
</td>
<td>
RunscConnect
</td>
<td>
ChangescConnectDB
</td>
</tr>
<tr>
<td>
<a href="https://github.com/SCA-IRCM/SingleCellSignalR">SingleCellSignalR</a>
</td>
<td>
1.2.0
</td>
<td>
human<br>mouse</br>
</td>
<td>
RunSCSR
</td>
<td>
ChangeSCSRDB
</td>
</tr>
<tr class="even-row">
<td>
<a href="https://github.com/sqjin/CellChat">CellChat</a>
</td>
<td>
1.4.0
</td>
<td>
human<br>mouse</br>zebrafish</br>
</td>
<td>
RunCellChat
</td>
<td>
ChangeCellChatDB
</td>
</tr>
<tr>
<td>
<a href="https://github.com/veltenlab/rnamagnet">RNAMagnet</a>
</td>
<td>
0.1.0
</td>
<td>
mouse
</td>
<td>
RunRNAMagnet
</td>
<td>
ChangeRNAMagnetDB
</td>
</tr>
<tr class="even-row">
<td>
<a href="https://gitlab.com/sysbiobig/scseqcomm">scSeqComm</a>
</td>
<td>
1.0.0
</td>
<td>
human<br>mouse</br>
</td>
<td>
RunscSeqComm
</td>
<td>
ChangescSeqCommDB
</td>
</tr>
<tr>
<td>
<a href="https://github.com/saeyslab/nichenetr">NicheNet</a>
</td>
<td>
1.1.0
</td>
<td>
human
</td>
<td>
RunNicheNet
</td>
<td>
ChangeNicheNetDB
</td>
</tr>
<tr class="even-row">
<td>
<a href="https://github.com/tanlabcode/CytoTalk">CytoTalk</a>
</td>
<td>
0.99.9
</td>
<td>
human<br>mouse</br>
</td>
<td>
RunCytoTalk
</td>
<td>
ChangeCytoTalkDB
</td>
</tr>
<tr>
<td>
<a href="https://github.com/SunXQlab/scMLnet2.0">scMLnet</a>
</td>
<td>
0.2.0
</td>
<td>
human
</td>
<td>
RunscMLnet
</td>
<td>
ChangescMLnetDB
</td>
</tr>
<tr class="even-row">
<td>
<a href="https://github.com/Elisseeff-Lab/domino">Domino</a>
</td>
<td>
0.1.1
</td>
<td>
human<br>mouse</br>
</td>
<td>
RunDomino
</td>
<td>
ChangeDominoDB
</td>
</tr>
</tbody>
</table>
</body>

**Note**: Some methods also support other species, but not yet
implemented in CCCbank R package.

LR prior databases
------------------

There are **13 built-in LR prior databases** in CCCbank R package:

-   CellPhoneDBLR
-   FANTOM5
-   ICELLNETDB
-   NATMIDB
-   iTALKDB
-   scConnectDB
-   CellChatDB
-   SingleCellSignalRDB
-   RNAMagnetDB
-   NicheNetDB
-   DominoDB
-   CellCallDB
-   CytoTalkDB
