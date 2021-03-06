#!/bin/bash
FOLDERS="
$MOLDB_PROJECTS/Clusters/Methods/c0/QM/F383118AD5D4ABB691FF0AA2A4B2E4A1
$MOLDB_PROJECTS/Clusters/Methods/c2/QM/24CA91D826028DEE38C49AC0301D155C
$MOLDB_PROJECTS/Clusters/Methods/c3/QM/D08D43EBB43687D0A0C5123EEAEBABF2
$MOLDB_PROJECTS/Clusters/Methods/c5/QM/AF807A6E060721995E9D081DD6FFB526
$MOLDB_PROJECTS/Clusters/Methods/c6/QM/FF83F789106751372FB75EC245AD8658
$MOLDB_PROJECTS/Clusters/Methods/c0c2/QM/12A71D6E3990A3367B970F137AFCE623
$MOLDB_PROJECTS/Clusters/Methods/c1c2/QM/B56F8D316913173260050B54EBFB4CB6
$MOLDB_PROJECTS/Clusters/Methods/c1c3/QM/7141054D110FF447F96483619434D693
$MOLDB_PROJECTS/Clusters/Methods/c2c3/QM/6DC471D6B2579402A2902E647C686977
$MOLDB_PROJECTS/Clusters/Methods/c1c4/QM/A8DFD71A34861DA6A26566819CF1C17E
$MOLDB_PROJECTS/Clusters/Methods/c2c4/QM/372BF74B8F68FAE75D48E10F9C36D84D
$MOLDB_PROJECTS/Clusters/Methods/c0c5/QM/FA318C3ED0E5C97E90DC90E2CC309735
$MOLDB_PROJECTS/Clusters/Methods/c1c5/QM/0C0A7FD47B6311B23A46219C7DEDA8F6
$MOLDB_PROJECTS/Clusters/Methods/c3c5/QM/022C36F7D5D978FE1659980C6BD82662
$MOLDB_PROJECTS/Clusters/Methods/c4c5/QM/9D0B0AF201E194BC311155D5472FCBE3
$MOLDB_PROJECTS/Clusters/Methods/c1c6/QM/EB3B84E0C1C84CECE43E8B379FE5E29F
$MOLDB_PROJECTS/Clusters/Methods/c2c6/QM/6E88CD2DA73DA17EE0035F901D7DC6FA
$MOLDB_PROJECTS/Clusters/Methods/c4c6/QM/6797D29F64488D308D192DAB4639C152
$MOLDB_PROJECTS/Clusters/Methods/c5c6/QM/373943C84F13B3926DA50C8D4349BFB7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3/QM/A2D6CE99FFE212243BFDFDDB0E12CDF1
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3/QM/CA651C56E72D20A70F416514B8FE3C7B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4/QM/6DB02CA9C61313AFD56BACF8717390AD
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4/QM/323EDE62F1AF8B04DFA257D91512D8DF
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4/QM/AC1A03946C843B23A42F8ED5BED0BF99
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4/QM/4B1378B7BFDEA4A56D6EE11FF52A2FE8
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5/QM/F3A1795854E3582F02B92D4E87BBEADF
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5/QM/BE3FFD8B198716FD929FA5FDE270A146
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5/QM/533C7FCF87AC19652DC1AD99398F6B16
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5/QM/7129BC3B8D430C89C3B2F3C240067323
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5/QM/E7BE2BB69DCEE9187374A0B2EABD9DC8
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5/QM/80F4CC4FE5395D7E40AC781536B4A108
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5/QM/5A19A515B639F4D9C8454DBCA89A903D
$MOLDB_PROJECTS/Clusters/Methods/c0c1c6/QM/701C42D0E6BC5910A172C45EA338455A
$MOLDB_PROJECTS/Clusters/Methods/c1c2c6/QM/3751392837C624A01421CB1D9274EB0E
$MOLDB_PROJECTS/Clusters/Methods/c0c3c6/QM/C1D15AB5B847CB34AD81DEA84150691F
$MOLDB_PROJECTS/Clusters/Methods/c2c3c6/QM/9E13FBA24E081F626819867E741E1780
$MOLDB_PROJECTS/Clusters/Methods/c0c4c6/QM/385FD44B42EE223DCE9D1C6D190D907D
$MOLDB_PROJECTS/Clusters/Methods/c2c4c6/QM/2A72EB62EB0B2B6DBCCABD1445BF319F
$MOLDB_PROJECTS/Clusters/Methods/c3c4c6/QM/3BDEC822260AF2045BD25972926094BF
$MOLDB_PROJECTS/Clusters/Methods/c1c5c6/QM/B9499E72A97C50780D0C131EC0A12C97
$MOLDB_PROJECTS/Clusters/Methods/c2c5c6/QM/8A86145AA8C28D028D93CCE876BF2B36
$MOLDB_PROJECTS/Clusters/Methods/c4c5c6/QM/44CD21C8ABA875F67F748BD42D444B73
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3/QM/E4D427B2D6E16A64D659C0517EBAB9A7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4/QM/BEB7708F10C22599F460B9A69ECF4A65
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4/QM/F279E6E3E63B1908B6CF2A50C0B81B00
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5/QM/6ADCA807AA480B6EA469048C4390A0E6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5/QM/9B167EBF1219802BF037C442736D3923
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5/QM/72910AA71C6A4546F3F6DA5161950CD4
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5/QM/4582CFE2A2BAB43A544EA474466BB6A0
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5/QM/D359C008CD61162F23B0416FAAFDA04D
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5/QM/C5CF34E812ADDA4C2A3E4D470B1DADA1
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5/QM/106F0E781A5106CD49DE0FE4A668233C
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c6/QM/91B8D542E0E75D2DA48C76D3EB422660
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c6/QM/5F15CAD9785366A4216C1C3DD2074444
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c6/QM/D412908A82E0454B181974DCB667B2EC
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c6/QM/248AAE30D768ED96466F3A6E0736A90F
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c6/QM/09CB3ACE1E9253B661196198CD580868
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c6/QM/4C08F1A727F1584C2EC332106BE08A99
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c6/QM/384EA1B2942E93E8C7647F8E2B2B71C4
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5c6/QM/80A3C5F4E6E2244C640D0999335E6A99
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5c6/QM/FACB9A099E987590291B7D036244E331
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5c6/QM/2E715FD99FABE25C12C603B59894ADF9
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5c6/QM/9E35534CF20DCB1A57BAC98835F05E3C
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5c6/QM/B91EF19DD542E19431E00B389FA0D2BA
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5c6/QM/DD9D1EF2B7B3FB15F62AEC6D6A06CF8F
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4/QM/31B14026162B9CDAA40C67FD69835644
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5/QM/0F2C0EED755F42D30C9D6BECF74C5A9A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5/QM/3C6CD093D438273080FA539C961B1CE8
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5/QM/47DF2CBE5A78350E704AC08FF6173B33
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c6/QM/92FD6991B4A1C9E65AD4AA1CDC56320F
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c6/QM/D0ABBBA7EA1C6D5B712279F4FEB32799
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c6/QM/267FAA87440CCB0931FE98B703168B50
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c6/QM/C7E8581FCBB866719EC108ED2C973894
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5c6/QM/7AF356E53AAABBA12B7A9F5CF0888EE7
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5c6/QM/D4346B4CD8AFF454E7F0B995EBAC2C67
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5c6/QM/3B4A3AA3B18457AC047AFAB3BFCB7E96
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5c6/QM/F25BBE4B64ED1309728C866DC95FEF6B
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5c6/QM/676FDB0CA628F6032335F319F3FCAC6B
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5c6/QM/C9063B8FDED5DF692C8929BBFE661FD8
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5/QM/2ABE44F0548B8424F2DFB86C8CCBC025
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c6/QM/718A6FB93BBAD984B4DD0FB17AF9AD55
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5c6/QM/367884A61CF444A2057115C080A48453
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5c6/QM/4BDECAB4BFCD7939BDC09AC9ABB64C31
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5c6/QM/C189BAB81DA948A691859579F97D19D3
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5c6/QM/D22F949A012A289DBDCBC8BDE4B51686
$MOLDB_PROJECTS/Clusters/Methods/c1/QM/328DA9A8E84519B52FB51C4543D08138
$MOLDB_PROJECTS/Clusters/Methods/c2/QM/38EC36FC7E07F1D918BF7F19905B307B
$MOLDB_PROJECTS/Clusters/Methods/c4/QM/D0AC7665AF09618B67443E687ED041D4
$MOLDB_PROJECTS/Clusters/Methods/c5/QM/950D094B9CB32F31FEB5DC063E8B559C
$MOLDB_PROJECTS/Clusters/Methods/c0c1/QM/D406F57285F0253D2EA58BFF155FB35F
$MOLDB_PROJECTS/Clusters/Methods/c0c2/QM/18A70890E20CB1B0D1454998C0971358
$MOLDB_PROJECTS/Clusters/Methods/c0c3/QM/B92A3A3A86BC1E3402836692A4B3E6B5
$MOLDB_PROJECTS/Clusters/Methods/c1c3/QM/20704FFA7A23CC0A2B2421E1F37045C5
$MOLDB_PROJECTS/Clusters/Methods/c0c4/QM/B10A0956247686D7227A3E0C7D849200
$MOLDB_PROJECTS/Clusters/Methods/c1c4/QM/22287793FE396FF2F3A8985E64564B3F
$MOLDB_PROJECTS/Clusters/Methods/c3c4/QM/5E6D4F83D28AC412C947ED53C4F628A4
$MOLDB_PROJECTS/Clusters/Methods/c0c5/QM/E22DFA48BDC7CF18C4BCE7E21E1C17C6
$MOLDB_PROJECTS/Clusters/Methods/c2c5/QM/7C5EF1E003CEE66420EAB96C372B354C
$MOLDB_PROJECTS/Clusters/Methods/c3c5/QM/747C012FEF8584423DC25F6E4E13B158
$MOLDB_PROJECTS/Clusters/Methods/c0c6/QM/9117329C15D136F7240450E79C135329
$MOLDB_PROJECTS/Clusters/Methods/c1c6/QM/4AEC6F6CE4DB059C5F69DCF33302AA27
$MOLDB_PROJECTS/Clusters/Methods/c3c6/QM/7A02B27F7DA618AF697D425806D19905
$MOLDB_PROJECTS/Clusters/Methods/c4c6/QM/3FBA8945F4CA62FA98752367012A53B3
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2/QM/55339D66485F2DD083903FC64BABC137
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3/QM/80B23EC42E85CA7454E0923C7A950CF3
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3/QM/A3212125CAB9576E11C63A280ABBD94B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4/QM/845A1CEC6F821D8BCE67037324A1B5A2
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4/QM/24FE8C7758C754EC72D40F75C8D54851
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4/QM/8A8D7F103B52A99C0E0519F240F19EAF
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4/QM/31F4FD8008EFBDC742E070FED41A7178
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5/QM/874DE80113969ECDBAEC49A06715210D
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5/QM/C274F3406307179C4ECC8726B36D2521
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5/QM/AD9DFA2FC46851ABFFB8BC94E21C9B64
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5/QM/BF59AC0A0E1662D8C70C0AA6B2A52D5C
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5/QM/664A23D2AF0A580FB2ACF5A5B94AC944
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5/QM/8017F4A651E51EBEEDE434A7C91E0CA0
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5/QM/4DCEB310A1BF9719698F212651183FB6
$MOLDB_PROJECTS/Clusters/Methods/c0c2c6/QM/C6948F752658EC07E9E11FF18BACB742
$MOLDB_PROJECTS/Clusters/Methods/c1c2c6/QM/6B1588D2D27D62906A62381606CF219E
$MOLDB_PROJECTS/Clusters/Methods/c1c3c6/QM/CCB99E703E96181F8D910205107497CD
$MOLDB_PROJECTS/Clusters/Methods/c2c3c6/QM/EEBF647E752F27A7129272B5E5184F82
$MOLDB_PROJECTS/Clusters/Methods/c1c4c6/QM/FAB5F9EB9C0F4B6526A0D9AB86B0DE20
$MOLDB_PROJECTS/Clusters/Methods/c2c4c6/QM/98DFFCB0FE43C546DB2DBC8AC3E60BF1
$MOLDB_PROJECTS/Clusters/Methods/c0c5c6/QM/F12E5A943FF9338D7AA508E9E0A4B052
$MOLDB_PROJECTS/Clusters/Methods/c1c5c6/QM/22FA00248027886282A7C1621FFD0B71
$MOLDB_PROJECTS/Clusters/Methods/c3c5c6/QM/6D012A0730EA363871903F7C70DD3A3E
$MOLDB_PROJECTS/Clusters/Methods/c4c5c6/QM/F0F6AB202AF75981C4E298FB4D0F2F7B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4/QM/578FD916DD3B95E628D9F9B8FF9DA642
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4/QM/DF1D269742BC5CF2CCFCAEBA759B92B2
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4/QM/59580B2C5BBBDDFA5494F6468477FA93
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5/QM/4095652D58D3B7D270455CF7371539F7
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5/QM/AF29F26ABBCDCE8B76EC95D36E1CAFBD
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5/QM/88F9A82CB9A80C56D98E027A8D78E3C7
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5/QM/81D94749FC3E86B139636A1566882005
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5/QM/E3A3B19B070659FADED297CEE3D62038
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5/QM/ADB87885A3AD57A800C4E6536E00AA45
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5/QM/03119406A32C79E80BAAB45E831AE646
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c6/QM/C1B128845502B6E797829FEF898A9344
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c6/QM/7BE7CF0DCB31EC0BDCC06FEBFED2FD96
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c6/QM/F78102033C5A9999B013F0EAD3CD4148
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c6/QM/8002195A93A15DD553ED76C3E0BBC68D
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c6/QM/96305B940F8906203B093B185E2A7716
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c6/QM/90DE1BAEAFEE58574803BC7EE8F97339
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5c6/QM/F741C107CD96B8DD00E1765E206E64E0
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5c6/QM/684FCE4AC7927C783DFC89A8EDB86A10
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5c6/QM/7BB5FE9B7300A0C2C652DFC0967E7DFC
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5c6/QM/98D1B2AE84B12B0C4267DA4D3D44A343
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5c6/QM/61B8EDDEF1ECEE748A7B586663251EA9
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5c6/QM/0A7EC941C195FEF477428DB3FA1AC9EB
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5c6/QM/E55A002BE0A2AE61970A2F765F33DDDE
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4/QM/C3FD2E4D5B796CB4C7B33906B86F73E4
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5/QM/B8AA10FFE70F4BFF594A742DB158DBB3
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5/QM/081F7B2B6477922BEE27DA2DF084023E
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5/QM/9CA91A951A268DBA31531EB2AC04DCC6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c6/QM/5E13EA81693DC69D1746A056BB4FE037
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c6/QM/F0D9E372FA8A520B73CDBCF88FEF3400
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c6/QM/9D3FD5B107DC932D52396DC92C282AA4
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5c6/QM/6411E50276A2493A72B4858D26D53175
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5c6/QM/86824422E48CFF8B0022C70D28E58F58
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5c6/QM/A895C906C80F9060648476EB634FF3B9
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5c6/QM/E18C77704236E670617E00C193BB6064
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5c6/QM/C4CD3961DFA899F835A2473BB1A68D5B
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5c6/QM/C80C0403856468F936709A2A9A536F27
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5c6/QM/7D4CEDCDD442C432D4856A6D427C5183
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5/QM/6C30F29FE0EC42FD9C96590D4CD88918
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5c6/QM/CB37419F2D9AE7F292FB8DBFFEA473BB
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5c6/QM/15F9C9F2236462A4D2F2FCE0C2290000
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5c6/QM/669EC5DDAEE7DCE8E95C99ABB9B4BD53
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5c6/QM/ADEFB34CAA67498AB1DDDF6B2E0CBF7A
$MOLDB_PROJECTS/Clusters/Methods/c0/QM/E314A3D6C7D70478DBD305FE44B93BC8
$MOLDB_PROJECTS/Clusters/Methods/c1/QM/182C8FD45EC9320042DEAD9F80E71E01
$MOLDB_PROJECTS/Clusters/Methods/c3/QM/272D60B474B2EA5E79618DB26E2AF292
$MOLDB_PROJECTS/Clusters/Methods/c4/QM/756B6E00FED19795E30283478C690762
$MOLDB_PROJECTS/Clusters/Methods/c6/QM/F0909560ED06A7C02DE3135E6B4FA2D1
$MOLDB_PROJECTS/Clusters/Methods/c0c1/QM/890A8CE3AF47B96E39D2F5B7A750CA4C
$MOLDB_PROJECTS/Clusters/Methods/c1c2/QM/7FEDB8C9146A5AE2146DFEF0BCE3B734
$MOLDB_PROJECTS/Clusters/Methods/c0c3/QM/3600815CB045D1E78951A5807B300BE5
$MOLDB_PROJECTS/Clusters/Methods/c2c3/QM/4B4E7C77DE4710B6A7CE3DBDC3B15317
$MOLDB_PROJECTS/Clusters/Methods/c0c4/QM/B7FC8E4D58971274C6881DA8586DD688
$MOLDB_PROJECTS/Clusters/Methods/c2c4/QM/4553D669EE1791AE85027FB89FAD5BDA
$MOLDB_PROJECTS/Clusters/Methods/c3c4/QM/F32619288C7966E841DBD5BE59BBE9ED
$MOLDB_PROJECTS/Clusters/Methods/c1c5/QM/5A404EF2BDBFE20920BCFAC3FECCB96F
$MOLDB_PROJECTS/Clusters/Methods/c2c5/QM/554B6D17E96BB75140BF86EEDF04E3B9
$MOLDB_PROJECTS/Clusters/Methods/c4c5/QM/8A3AF3C86E92A6BD9CFCDF0F04E65805
$MOLDB_PROJECTS/Clusters/Methods/c0c6/QM/3920F77327648544F0389FC4C24F50A1
$MOLDB_PROJECTS/Clusters/Methods/c2c6/QM/3904B25B6309BED3618586CF75292BD3
$MOLDB_PROJECTS/Clusters/Methods/c3c6/QM/C59293D62FCCB069967FA543C730D747
$MOLDB_PROJECTS/Clusters/Methods/c5c6/QM/DA79D1B5C77CCC888CA8D48ACB3BDE0B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2/QM/C1CB2D3246F997263C45FDC407B64B22
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3/QM/2F164C18AC3478DAB8E30D2AE8455681
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3/QM/80FC34DDADB1C9F26D95B78C942B8F50
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4/QM/756F6E4D37019B467DA47FA1947CC106
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4/QM/0967BB072C7DFA4A569BC91F4A865A79
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4/QM/9093B42E81A53B6DF3641B4230D7CB02
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4/QM/8DCFAE0ED87614DB13A2CE3E938B0CB4
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5/QM/D18E1A9201CC32F60483A3FD3A09CB02
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5/QM/5F3FDC941C351EA880C552F52E4074C5
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5/QM/007BCB880D485001A2B51B9A27667503
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5/QM/B7ED1DE69A3A0BA587537F85440B89B9
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5/QM/58CB9574A250C5561109DEDAD40695F7
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5/QM/ED30161BA4A908E7330166AD7619BF59
$MOLDB_PROJECTS/Clusters/Methods/c0c1c6/QM/77A26440148F5D819E1285D035CB9EFD
$MOLDB_PROJECTS/Clusters/Methods/c0c2c6/QM/2878868448D16AD8D1306B3BA3E8FF89
$MOLDB_PROJECTS/Clusters/Methods/c0c3c6/QM/55BBE3EDABD8DC54BCB8E0284F3E1F81
$MOLDB_PROJECTS/Clusters/Methods/c1c3c6/QM/E9D61B7DBDEAC38B4AB588F180ACDEA8
$MOLDB_PROJECTS/Clusters/Methods/c0c4c6/QM/3EA37656464D413917B533C63FCDDCA5
$MOLDB_PROJECTS/Clusters/Methods/c1c4c6/QM/F94B5DC8FE1713B5617252BCE79B6719
$MOLDB_PROJECTS/Clusters/Methods/c3c4c6/QM/FA9A92B26952C1D2D0BBE9AA7C04410D
$MOLDB_PROJECTS/Clusters/Methods/c0c5c6/QM/8A3E08941AEFF132E8DC051F9DB134B8
$MOLDB_PROJECTS/Clusters/Methods/c2c5c6/QM/B49356D5D8EF2EDF40E025E7FDEF0D15
$MOLDB_PROJECTS/Clusters/Methods/c3c5c6/QM/16C9F07B29B3EA779662F7CFC820456A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3/QM/6C2A404960857646B74800348276EB8A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4/QM/A6F86DB29BB4BEDEB04CC6C33C612FF4
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4/QM/8DEDD80BA08C881ED6414DFD140A2EA1
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4/QM/FAC49707FD1E1DB356FD37DD2154A895
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5/QM/30B45F1E3EEA54611E8C6FE0FC4A3D98
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5/QM/E3C36F7C7A191C0E13B01D908EFE63B5
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5/QM/A69EA5D1DF650A3F70E47952EC78F48E
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5/QM/CBDAFBEE4646891A94E893AA38A78B46
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5/QM/DA75B0C3DF8744D1B109FA592F91D59B
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5/QM/FBE51B47775B6BA8255C80F74F548BE6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c6/QM/A36DCC2C5612DED678974419EF7AFAC7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c6/QM/DCDAE687F913A3EE51872BF30D7760B5
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c6/QM/995F66004EE17B9FD7F844DD00E1300A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c6/QM/CE701BC5321BEBA8CDCDD369310ED458
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c6/QM/7CB4D106960F7B32AF7CAFF76967E1B3
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c6/QM/3BC8F990DC72C9614EEFBA63859400C0
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c6/QM/CA317E31C30B62E25F015D6E563E25D3
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5c6/QM/114252BD1B50AF09C9F7AB11500285AB
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5c6/QM/55F8228AB7AFC6CB488540CEEC745361
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5c6/QM/CA2FBDF58734916EF6D480FDE23742DC
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5c6/QM/13DD52E410028250A346378A22490034
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5c6/QM/B74BF4CC26D72637E607C8DDAC650126
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5c6/QM/B47391DBDC4D6D517C98E9985FCAAA6F
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5c6/QM/FCA1D0499BC8BDF764B4B7A9187F7595
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5/QM/7293C2EFD4CF49D13CA13B72101EAF8A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5/QM/B8F046A9BD728327E4CAAF1D9E9B065D
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5/QM/EE552F44A59FD9E63F2FC5B48626788E
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5/QM/9DD9606D9431AD7530D338D7F555D843
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c6/QM/7CC436D6A1DD63DDDCC71EF5F5B3602B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c6/QM/57313BAFA6DA92B946F10129BD7721D3
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c6/QM/8A9B12AE306070DFFA49814D6809E569
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5c6/QM/7CEF92548A23FBCC185F5814A8CD0B61
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5c6/QM/76A90821CC9C65C2DBFF3CCD6FBF7061
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5c6/QM/E6C7AFC4C8662F2FA288BC8BB97CE809
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5c6/QM/70C59189B8D4860B1ADE5B7576D2EC0B
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5c6/QM/F4D10216C27F574D5C30FA4C853E0759
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5c6/QM/AD0B715EB626C636299BA76038AB4664
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5c6/QM/0F17BDF8B1EFF66467C87AF795A8A856
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c6/QM/0F5FD50C6634F6CB595C9D24FD0912AD
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5c6/QM/B102416FDDE2AB9D065CC0C0CE89E3AB
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5c6/QM/64D70EDF3DBCEC9F15C8097135F64AE9
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5c6/QM/751651DCE3E15B925B8BA31737E6A679
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5c6/QM/10565F8F6E730B80516CAECEE7002C49
$MOLDB_PROJECTS/Clusters/Methods/c0/QM/3D82DEA8836B57C40EB2A9B1EDCE3789
$MOLDB_PROJECTS/Clusters/Methods/c2/QM/A361C881F127DC1F7FC5BE3FACEF2F90
$MOLDB_PROJECTS/Clusters/Methods/c3/QM/B0075572E28D19A97340F2F3C36593CF
$MOLDB_PROJECTS/Clusters/Methods/c5/QM/409E5D0BABBC78967E44D1076507A9C0
$MOLDB_PROJECTS/Clusters/Methods/c6/QM/725BD7E579A6D24CCD1A48B0E32E620E
$MOLDB_PROJECTS/Clusters/Methods/c0c2/QM/329DC1623B4BC8F29737669E7F482C01
$MOLDB_PROJECTS/Clusters/Methods/c1c2/QM/C8FE731BD1894D9FF89A2C9150AB6C93
$MOLDB_PROJECTS/Clusters/Methods/c1c3/QM/3E29A7382A497FBF82EBBAA451B95337
$MOLDB_PROJECTS/Clusters/Methods/c2c3/QM/4B2B19907D7E200BD5244894E221EAA8
$MOLDB_PROJECTS/Clusters/Methods/c1c4/QM/2663C30A80CB5FC1FB2FF125B42BC705
$MOLDB_PROJECTS/Clusters/Methods/c2c4/QM/AB22072B1AB709A5A2C96B1628F25D1E
$MOLDB_PROJECTS/Clusters/Methods/c0c5/QM/CEF296BB6E1FB3B9574527121836A9F5
$MOLDB_PROJECTS/Clusters/Methods/c1c5/QM/DFA5898C9FB94A18B400EFD354DF736F
$MOLDB_PROJECTS/Clusters/Methods/c3c5/QM/C1877B63779231A5C0AAFB68B1616665
$MOLDB_PROJECTS/Clusters/Methods/c4c5/QM/67304781C944C04F555225D19BB063D6
$MOLDB_PROJECTS/Clusters/Methods/c1c6/QM/52B1BACDAB5A155ACC7E33631771CC55
$MOLDB_PROJECTS/Clusters/Methods/c2c6/QM/121754704D42C377C4E25E417102F7E5
$MOLDB_PROJECTS/Clusters/Methods/c4c6/QM/938B845DDB29D545F570689B4ECD6FB6
$MOLDB_PROJECTS/Clusters/Methods/c5c6/QM/BE0FB937B27FC9612A41D2C8270B2DE3
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3/QM/44AEEFB5A1C4A8E07E64347691B6E9E3
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3/QM/325254D88C5393A5381E72BA46D8AC46
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4/QM/306B3A77B24D85F53E5BB58EA2C206F8
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4/QM/2709C70558A4949AAB6DCF24763B317D
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4/QM/88A97EBF2E06F6EEAFA83E6215697287
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4/QM/B8DF09CD25BF7A8942AA155A536D421E
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5/QM/BD643B1CDED41A0A32640E8673DFC9C7
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5/QM/F02F0E2FF8A2A7E4CE8289B98547CEA3
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5/QM/4AF9CEC025DF6A90E7B81B354B199F00
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5/QM/24A32E6AA33D72CF9C69D2633526FF97
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5/QM/D809457F6033106A2AD72185A8A173CE
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5/QM/95231AB0541362052203CE31B9EC152F
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5/QM/9A8B33A00E36CFAA4BE7BB8FDA78A3BF
$MOLDB_PROJECTS/Clusters/Methods/c0c1c6/QM/CC1D207CFA63CA3A14BAC8F6B3CE6F74
$MOLDB_PROJECTS/Clusters/Methods/c1c2c6/QM/A971049DCC21AE3907F1222E0B17905E
$MOLDB_PROJECTS/Clusters/Methods/c0c3c6/QM/E538D6F65701723674D16DA3409C849D
$MOLDB_PROJECTS/Clusters/Methods/c2c3c6/QM/B243FCDFB268B42EF6AB6D0C979F2F74
$MOLDB_PROJECTS/Clusters/Methods/c0c4c6/QM/0249AF065913ACC2EFD5C9C20F3266A8
$MOLDB_PROJECTS/Clusters/Methods/c2c4c6/QM/B889D8759EF6079EFA113A6A51531922
$MOLDB_PROJECTS/Clusters/Methods/c3c4c6/QM/E4CF900764FD8440252DA2BF7647F7B0
$MOLDB_PROJECTS/Clusters/Methods/c1c5c6/QM/065360397F45C5A80C93DA37BE48672A
$MOLDB_PROJECTS/Clusters/Methods/c2c5c6/QM/C29025E9EC4E3A412923E8B01B0BB8C0
$MOLDB_PROJECTS/Clusters/Methods/c4c5c6/QM/BFE24E234494E0E8F4840EF58A2D3F15
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3/QM/C401FD309506A13167A92549EE765D4E
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4/QM/8E145821C8DF7F5CB747A3FA391F4F77
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4/QM/3F1126F7427A2D81A93ECA98D16C5FE0
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5/QM/9F556DB227CF5E041A2F30B3738B7CE2
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5/QM/F0F57DD080D2580451C04F377F8A1FF8
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5/QM/233988556DD46827783DBAC80762EB8D
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5/QM/BABADCA65C185C19A7E426F354C94CF3
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5/QM/0461136D6DB15CABE3151BFC9D55BEC0
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5/QM/59A3A92DCD091CA4A79709F6B138CE6E
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5/QM/AB8A3D4867A6677C1F186DEDE95EC5B0
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c6/QM/61871A9880EBA2996778E6C134FDD23B
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c6/QM/0B7D874B818646236DCCC3861E637921
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c6/QM/CAD7A84A4E18F5D12248ADEFA8630C06
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c6/QM/E45C88385D52C32B2062605B4CA26DF5
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c6/QM/8D1A72F0247FA55AFC899126941F2803
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c6/QM/BEB1169518848D9FC87F30BE0DFA43B5
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c6/QM/BBEB69EC0FBDC177C4B8FE4A446262BA
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5c6/QM/FFEC14E65FCA07582D6CC012DF87B763
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5c6/QM/00B906F4BD8D5DE24FF2A51E8675702D
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5c6/QM/B9F4625AD66F7A5FABB8056FEFC7BF40
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5c6/QM/354C1C8026566143F44482C5E5688F66
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5c6/QM/7DEEF7C46A77256F01F64E707282BE56
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5c6/QM/7F6DB2D7ABB67F89F2E6338B667B326E
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4/QM/8EB4DA03485514BAACB008659B3B2909
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5/QM/D0965B31C2A7F9A6F0FF8BC74C8A907F
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5/QM/91567C9C6C19ECB2DE68F9D65D6C2C16
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5/QM/48D4D08A7D5C71E999CEDE9B98C3827A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c6/QM/8939D66F95EB7845A5810BC6E93DAB97
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c6/QM/501C9FA7247F6C4CFE577FEA6568385F
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c6/QM/A5A08847D7427F1C1F3F30C767C7D2A2
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c6/QM/6A2925AA378CD1335D4EF74650D9B5A0
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5c6/QM/2C4F75FFCB83B4EFB98A6990CB79D4E2
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5c6/QM/D7CABE499BE3343F96F1A1B4345A7A70
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5c6/QM/70879D814665CCD35AD18BAAD5659A37
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5c6/QM/A9D3CACEEB7C64806FE083A48BEDB232
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5c6/QM/0DA5B11A66D957631E738CCD742242E9
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5c6/QM/8FEA71469933909747EFF2B8B0F476BB
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5/QM/2FBCE5FCBAD7836F8ADF5F438A9DBDF6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c6/QM/5E3746709F0DB93D83DF93698BEC2C7B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5c6/QM/C5BD0277D5AB6A75E437EF84DC3A222A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5c6/QM/B210E5E2C5A774A1AFA2C4607D440127
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5c6/QM/95A90855826D2F1B77DF9BD6E9F49BA7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5c6/QM/7776DAA02A71FA2121558775BB144101
"

N=$(echo $FOLDERS | wc -w )
count=0
P=$(pwd)
for FOLDER in $FOLDERS
do
	echo Molecule $count of $N "(" $((count*100/N)) percent done ")"
	cd $FOLDER
	. ./PARAMETERS
	. ./DEPENDENCIES
	./run
	count=$((count+1))
done
