#!/bin/bash
FOLDERS="
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323/QM/593A93F421930B5D1BE2D07873188992
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323/QM/4FB73A43F48549E3250DE9E0704BEB19
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106/QM/502339B8EF2871DEF31D88A5C3BA8A77
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106/QM/7B74E8980FBFBBEEE1550DC17C6C680E
$MOLDB_PROJECTS/Clusters/Methods/w1w4309w124/QM/CFEBA35D19F808A3FD4AE0623C5F8248
$MOLDB_PROJECTS/Clusters/Methods/w1w4309w124/QM/39FCFB8023EA62B3A43E40C82049D938
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w3205/QM/BC77B7B7E1D3CD4C68A4831EA8541203
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w3205/QM/4286C7CB0907D38DA1211CB8360FF6B6
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w5014/QM/C2716477518D91CE24F08162448FB65F
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w5014/QM/C03907A19134D4FC34AE5FA09C7B679A
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w2293/QM/3FAF9048982D4F7FF57A633CE02DFD7A
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w2293/QM/ED63B3D9A591726443FE3ACB461DEA74
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w3106/QM/80BEFB2ABD40ECD2B431CAF4EB3B91B3
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w3106/QM/2641915FECFCD1F52586A8856C2C25B2
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w4645/QM/EA3B1BA83A2497721CDD832579B6CE9D
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w4645/QM/52DA27F69C96F7314D1341CEBF794AF2
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w5014/QM/1F3EDDE990983B0D0DEA34B4633BE39E
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w5014/QM/D22A3489D6628A3FCA3D2F4A84748D86
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w520/QM/99844BBE3E6460BD7E052605A41C28AB
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w520/QM/3B7C5E775BBBFF2555855FD41955D4D5
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w2626/QM/9F38BABFB14CC7CC9DDAB65B965E7C13
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w2626/QM/51DD77652F17EAD354E5C3F0CBB447A8
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w5068/QM/027A4F991F8B218CA54FF6B7CD31A26E
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w5068/QM/F9B0BE495C4CA75FC75BCE2475A51E4C
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106/QM/1BC1B67968D7AC8FE683E3D66B7DA1EB
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106/QM/FE64D5F50D13068FD6204622107FC08F
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w2323w4309/QM/AB16C815CD297542DDA47156F7D0BBD8
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w2323w4309/QM/C89C2C0651B53B70AFD8ECB75AF55842
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3205w5539/QM/8D3A45C221DB1975751C9326CD0087C3
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3205w5539/QM/EC456520A15C316F29F5553DFDDD97C7
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w5539/QM/A16DFE2A97A76ED6941F0C0EF7A699EF
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w5539/QM/5578567E4E0507E8D8425C8BCB3ACF04
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w2323w6538/QM/7E9042700CA094B63CE795F9189B2BD4
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w2323w6538/QM/8A2C88DF1990BB1667E6985BE99A028D
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106w6538/QM/EBDD9DDE10935BC5E7B9747114220303
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106w6538/QM/EC4CB826DF49B0FFEE75BB3EAEF94FF1
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w4645w6538/QM/98EA06098D4B3EF098E7A459998FBA71
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w4645w6538/QM/7F3E8B3A89D2A073C3199CAA0F13C0AA
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w6538/QM/4D51392BE914F2A4AC66174EECDCC3C8
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w6538/QM/5F4B257A52F97820A753D42F821E1A55
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w2323w6550/QM/2591544743E29BC49BA6F709937512EB
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w2323w6550/QM/0A21D4507F8A702D32D023C939138C82
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w2626w6550/QM/5940840DE04A2D8F1276C71D64DAA332
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w2626w6550/QM/AFE9D62BCEE640F5326E3CC38D2B7AA1
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5068w6550/QM/82B06D5892B36BF70F5E225D9B2EB025
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5068w6550/QM/8B7510D08FD203EE4F8E174B66BDAC6F
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w3106w4309/QM/023B576E39B3BED4998E84E68B88EEB7
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w3106w4309/QM/BA60A16E9C64620A3A15F96570437AD0
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w3205w5539/QM/77FB5241AA3A2CA195C5E7C0C22FBB4A
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w3205w5539/QM/3DC007597EE4400A1C79B5925B1BFE6D
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w5539/QM/D1E2D902D39C35332D01A5135E8228BF
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w5539/QM/39A8EADFA1DBF03E9957B84BA9C2CE3A
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w3106w6538/QM/B6C9A281F9D816D58170B5A8576964D0
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w3106w6538/QM/0C3E798595BB1847D260367BCD09AA8B
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w6538/QM/450270A5A40503BE569E7CDE8436143F
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w6538/QM/D1038B7A5393DBB3BE49518976B3A22B
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w4645w6538/QM/0D02D6298038D1E905C786B04454FA36
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w4645w6538/QM/0ADEB933F0B10A6F6E7D5F09CC37FFFD
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w6538/QM/207EF0B74A361B791D335BB712E60D19
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w6538/QM/76AEDCC59AC4E23D47B3C94ED991C933
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w3106w6550/QM/7BB2C22E89AA2288247D4F2EA0AB9B3E
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w3106w6550/QM/7B5133DF07A87A029BE2D3BFA5602E97
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2626w3106w6550/QM/2D6EABFC76733F5BF0DC823218048F64
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2626w3106w6550/QM/5965FC07FF921AD0F9F69EA44CE723E9
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5068w6550/QM/A1D42CE57199383BAFBA3CB91C955B3F
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5068w6550/QM/F20FF7C38170B99FF43542186DDE46B5
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3205w4309w5539/QM/737892E5532FD4DA1D8B784BBB7D8337
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3205w4309w5539/QM/387FE9E08AC86AD8A3F8FBBD9FBA89BD
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w5539/QM/A2D1B9C271958D84A3CAE1121D2FB441
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w5539/QM/C879BEC48E05C69D6DE2C08F87BFBCFF
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2293w4309w6538/QM/A99D4CB3F78A7D504402D386622DA9E1
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2293w4309w6538/QM/33DAEC450C7EA8655DE81797DE4A07EF
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3106w4309w6538/QM/654E9FDC37D75B03E9176D3BD45CC5AB
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3106w4309w6538/QM/B3490EE89E2B8FE4E35AF1D55E8FB1EF
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w4645w6538/QM/379C38D5255A015F1754B5C31F07165D
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w4645w6538/QM/66694B74849BD892238F8A431F4280C2
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w6538/QM/A6C7A50B52CBBFCF9B757D146F5D7515
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w6538/QM/E6167A39AE373482DAB861432F7E511D
$MOLDB_PROJECTS/Clusters/Methods/w1w124w520w4309w6550/QM/A27CB853C6D832E76BA15F56A925F5C2
$MOLDB_PROJECTS/Clusters/Methods/w1w124w520w4309w6550/QM/EE9FC04349A3D6D6107A1D4327698484
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2626w4309w6550/QM/178A6629FCD5EE2F328D0602E3621286
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2626w4309w6550/QM/8D8C58BA1A636399CE6DBE29A5F037D6
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5068w6550/QM/5BF4A44495E26E2ED9C84A4AA8227E0C
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5068w6550/QM/C052228FE7DF6C554F300EF335FA9C01
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539/QM/E3E7157F7EAED2286C6D36D4708AB45D
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539/QM/B290A52CF6E339076003E7475BD633E3
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3205w5539w6538/QM/8029102E9EAC410251A2A580D445E57E
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3205w5539w6538/QM/9A3BA174A1FADBB3C81D9565601BE2AD
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w3205w5539w6538/QM/DA5EF470360CE3B41EE74CE698481EC8
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w3205w5539w6538/QM/22E3B4AD556AD516BF6156BA390CD1E3
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w4645w5539w6538/QM/C9C043A2B513A5777DE6156F25A737E9
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w4645w5539w6538/QM/51F51B5D1426BFBDEEAA69D272114E19
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539w6538/QM/C102E6ED96CAE5A0722336AF60E5C883
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539w6538/QM/57CF9A82D90025B37E2831D1BDC5D40E
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3205w5539w6550/QM/FA758F8B677494DAAED6323E4C1B5066
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3205w5539w6550/QM/21E58D4F37F95AD7CB0491A4F02F27E1
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3205w5539w6550/QM/F1F23A293EB4821C414602D1900A76AB
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3205w5539w6550/QM/9A3F5F2E0F1AA93D7B24E6993EEB3715
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5068w5539w6550/QM/25D98DC011257999790E55A86CE208DA
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5068w5539w6550/QM/F3F00C46B30734A74043C0E2B3BFCF13
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w5539w6538/QM/67897CA4D76BF3FDFDDA44E350A9AA77
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w5539w6538/QM/9B300FE7172C138E5454986C396B503C
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w5539w6538/QM/A0CE919D3FCFFFC483059A69E1F651E7
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w5539w6538/QM/7AC8E5382E718ED29F813871CDFAB75F
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w5539w6538/QM/0C65D1CA0A4B4659C6DBD54B29826717
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w5539w6538/QM/1857D89E2D0D32402C03D23CA9BB53A8
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5539w6538/QM/474A2F57097EEBD063605BE16046BFC4
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5539w6538/QM/0A3A20AEABC098EE4D26616FB7440FDD
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w5539w6550/QM/D8002C9EE49681C0D2CADC73CD9B0960
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w5539w6550/QM/A17AC1D98F1438098D4B861AC4A5F060
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w5539w6550/QM/164EDAA50DFD290340F9FD0F9EA83CA7
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w5539w6550/QM/E6EDD55AFB1E3FD20908747C3B8CC687
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w5539w6550/QM/C27546C466B8C7822F4C37BEC64EF1DB
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w5539w6550/QM/F32F7A73FCE27F9B12CC87F287C73F4E
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3106w6538/QM/24B35207E95A5576DCF2CA2E3B71E2F2
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3106w6538/QM/7943BC81840995D664FD687C2E36DAF7
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w4645w6538/QM/7698FB79B81F91221D0F9CC02EED2131
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w4645w6538/QM/911882E5D22735FABC57436B27102563
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w6538/QM/903E835198579B11187D70DD411FF012
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w6538/QM/C8304656542D6C192885EC806BC7E8A8
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2293w6538w6550/QM/B7C6E73262653E7B2E4C125EE7DA7041
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2293w6538w6550/QM/B38B1331508D9F23DBEAEB06AE0C0EDB
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w2626w6538w6550/QM/44A7136F5F77172694DAB85BD058EAB8
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w2626w6538w6550/QM/2CB6F1F5D436EFF79AB007ECC5082502
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5068w6538w6550/QM/20DA4F83A36821A7D3F2BD6B9D87C2A7
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5068w6538w6550/QM/15BE9DFE256818256C5EF1D38A945049
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w4645w6538/QM/59399A586C7022CEDDF72E4D9982CE45
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w4645w6538/QM/5C44419DCAF3E9D5310304BF25C7E5B7
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w6538/QM/F289548536C84448C091424333960776
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w6538/QM/5E0103C9343431872F3B701F5EA5B82F
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3106w6538w6550/QM/72868B17179C8EF2B6B9FC6574C58BB5
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3106w6538w6550/QM/69BCFD614C906BCC01841799BC239729
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3106w6538w6550/QM/1B263A1226D78917B7F5819E01004F1C
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3106w6538w6550/QM/FCF00D4E55E5B02AF475520BF8E14614
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5068w6538w6550/QM/5C459E90EA2C1F8B60FD5300EEE74DF3
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5068w6538w6550/QM/EA671715EF1A69F5C9C23484F82F5BB2
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w6538/QM/09E4577886D50D4A9BD918229F18B9F0
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w6538/QM/663673951296AFA3CCE7D2B199A42F86
$MOLDB_PROJECTS/Clusters/Methods/w1w520w4645w6538w6550/QM/A119E7599AF330B96340FB7AC527B04A
$MOLDB_PROJECTS/Clusters/Methods/w1w520w4645w6538w6550/QM/CC19EAE34B576FAFFED69CC0F34BD032
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w4645w6538w6550/QM/7198B2BCF205A04D87F8784257DD8A88
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w4645w6538w6550/QM/2C637A3E461336384901D662916F7695
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5068w6538w6550/QM/C40A930DB749DB2FB6813F2070391336
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5068w6538w6550/QM/56FA6FCCD05C507B9187D6CE93851EC0
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w6538w6550/QM/4842DF22F9E1AFECF6197C0D02506B88
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w6538w6550/QM/EAC6BF18379754908A25CC52B9AAA4DF
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w6538w6550/QM/71017FCB171A5461CF4286822AD67F59
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w6538w6550/QM/83CA15F32F6073A3BCB5449F28D6FFC6
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w6538w6550/QM/13BF68829A0C0F7AF27E7EC7588F2FB9
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w6538w6550/QM/DE8E54C70CBB76F00BD156EF42F90DC8
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2626w6550/QM/C6A81DD20E09B8D48BD6660824D224A8
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2626w6550/QM/FF89D92517ABFB4AB4D5B7EA02285AA6
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5068w6550/QM/A860AE1B6D85056D53C9BD2A3AF67184
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5068w6550/QM/0F00C06F5423C6EC0BB6EE852C1126BD
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5068w6550/QM/84A1517F14BCAEE612FDF478CA3708B0
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5068w6550/QM/881CCFD57B3BBCB118DB141C6A0E9D38
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323/QM/8E83B286F99511AA7D5B0B1F3A25AC15
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323/QM/6D484AB369ED81FE5D5AA56A94E9D9A8
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106/QM/694D30F8329298E51AE6BDB4BA6243CC
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106/QM/D107683BF40618EA886A5600266A67D7
$MOLDB_PROJECTS/Clusters/Methods/w1w4309w124/QM/CE26D9538A120C12E8AD4EC071AF510A
$MOLDB_PROJECTS/Clusters/Methods/w1w4309w124/QM/067C1C007C6703C28B39FA6DE0E12234
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w3205/QM/050FD6DF69A7AFE51E9AA0BCAA30B5A1
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w3205/QM/763A477329A0ADF698EEC948D083157C
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w5014/QM/08B930D30A50D27234732D96FBCAE67B
$MOLDB_PROJECTS/Clusters/Methods/w1w5539w5014/QM/312DF8CA9A0F9D5CE6E09F23B95CBD21
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w2293/QM/F0CC9D8D434C741C9C20C7EB960EADD7
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w2293/QM/3D4099489BF64F76DF520EF1315281F9
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w3106/QM/8B0FCBC2BEF533C10B68A3D5F1C78353
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w3106/QM/440C871D2262FB5C204B09BDA2E0AC90
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w4645/QM/13A186534A060B93DE8B58A9981DF2C6
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w4645/QM/DF9BCC61E7D9752B48F1FA3C5057F4EC
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w5014/QM/4AEE192C8B7FCFFFA6CBE2E502808317
$MOLDB_PROJECTS/Clusters/Methods/w1w6538w5014/QM/CE939B2E44BFD05225E799EB507A52FA
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w520/QM/FE3D16CB8065FFD5F3954C3E3D13FD44
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w520/QM/E067F0B74EBD924FD0ECD9AA6EB99504
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w2626/QM/C374EA2182EDC8B6DF67217E86EFFC58
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w2626/QM/13B6FAB52FC84732E3832F51DCA66E7A
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w5068/QM/2B75E1D3FDC00B6C1E9143757093E762
$MOLDB_PROJECTS/Clusters/Methods/w1w6550w5068/QM/7D0888110D1360ADFF63B598650CB75A
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106/QM/1FEE08D60F7F0CF24F5BE1D6DF8FEC88
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106/QM/52EA5A3A55E68AB6735FA2E71E158062
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w2323w4309/QM/E34F8DB7D7725E6BF1E8521BA80294D7
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w2323w4309/QM/EF73C5ABC384FC2BD710BF13259CE71F
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3205w5539/QM/7870189595BBB3ED4ECB697D3ED8B982
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3205w5539/QM/E71F70103B9ABDD757CF77141E09CBDF
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w5539/QM/0639AD248CA02F2BC768B0D774E3E1A9
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w5539/QM/C230D6FD08A0DBE52D2F70009820EDFC
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w2323w6538/QM/4351AAE63A5CF8A17AF0A5C2F18CC357
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w2323w6538/QM/D0FF5662E520A31ACBACE906F0039C7C
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106w6538/QM/472546A28B782F0F66360FE71BE1672E
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w3106w6538/QM/65036601CFEE74013F3AE807371A67E3
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w4645w6538/QM/54C8EA96FF9BA58DF4290C7EA7C72B78
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w4645w6538/QM/A7C7A63C59AB4008600525F1568729D3
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w6538/QM/0A9CA012EB4019B655CE8A4B361F0E25
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5014w6538/QM/C7E7F6341D19A4B52E043C46665853EB
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w2323w6550/QM/A6F4A0034E16967B16A64D2C02E5E2FE
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w2323w6550/QM/4621D1D4944692AFC4A1606D26E73767
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w2626w6550/QM/253415B6BD0E45511FAE99FA0FF41EBB
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w2626w6550/QM/2048EFAA5F3343B81E425C36DD810473
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5068w6550/QM/DAF4F5E206F176AF83FF92582A640EC6
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2323w5068w6550/QM/3857DE4B6BB248BE7CE5756F5480EE14
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w3106w4309/QM/B1BD1AA5699A6DFF370E7280AAAAA7B5
$MOLDB_PROJECTS/Clusters/Methods/w1w124w1183w3106w4309/QM/DDA83D22E081FB4E1A22CFAEA6A48FE2
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w3205w5539/QM/C748D9CB9BB526B770F9D1F3B4FD46A2
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w3205w5539/QM/C45E29EE3DB2BECE1A7BBB2C18E0BF47
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w5539/QM/527111D6C3D6FF4AA38F28BB004EFFD9
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w5539/QM/20E2BB94AE7589638442CAB3F92CC13A
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w3106w6538/QM/99137C25948C49344654B9D319057000
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2293w3106w6538/QM/1E33586A1EF80386F9D3AFFAF9423413
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w6538/QM/8F9EF1C45D11C6E30A9F0AA8510DB7A4
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w6538/QM/142CA849C14234EA18F15EA544EA54D4
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w4645w6538/QM/33BA94ECAE7F54CAD882B4A7BD75270F
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w4645w6538/QM/F8568CF2509BF2AE73E2585A222912A6
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w6538/QM/401E98752DB1149CAA019000D1079156
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5014w6538/QM/0982A630FE8870006984DB4A80073692
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w3106w6550/QM/0FFC02AAC9AF5A11C9D73130A86624E6
$MOLDB_PROJECTS/Clusters/Methods/w1w520w1183w3106w6550/QM/96007B7CA88E56362BFD2A671AB36A9F
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2626w3106w6550/QM/6D20775664B3219F62E7273B2282C742
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w2626w3106w6550/QM/1D2CB9ECC93EF9F341DCB8AC019291AC
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5068w6550/QM/24714FF41E120A63F17F8155C23AC7D3
$MOLDB_PROJECTS/Clusters/Methods/w1w1183w3106w5068w6550/QM/9F156EF080A8B10D8A7B68DE156EBF4D
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3205w4309w5539/QM/9BF471FBDB08EBD2E0D434E55BDDC781
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3205w4309w5539/QM/234B57EF7AD33923C81EA5EBE02560FB
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w5539/QM/AE17A763FFB94B4EDF792DC993DDED5E
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w5539/QM/24C9F4CF9BA2CA2CA77F969729978394
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2293w4309w6538/QM/DEF253B893EE1B5AEF5275DB8F477543
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2293w4309w6538/QM/FCB24504521C5D618A617CE609F2A28F
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3106w4309w6538/QM/3749B3C0FC2EB72C297032A183ADEECB
$MOLDB_PROJECTS/Clusters/Methods/w1w124w3106w4309w6538/QM/D945ED04EF002B31A18A6AACB09C6939
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w4645w6538/QM/C075D0038F5C92F48F00A3080783BCA7
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w4645w6538/QM/96F9312AE99D2C37D4D4180AE9F43170
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w6538/QM/739AE642191B17BEDD300D2719223C20
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5014w6538/QM/B560F689DFA7CA816ED9ADA9DE7E06EB
$MOLDB_PROJECTS/Clusters/Methods/w1w124w520w4309w6550/QM/5FA0845CEB0597860D652EC2A0E16B5D
$MOLDB_PROJECTS/Clusters/Methods/w1w124w520w4309w6550/QM/7A6A2195775BCA6D4C6FCFFEB05772CC
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2626w4309w6550/QM/F2DFE7AC69F8796BDD689ED15167008E
$MOLDB_PROJECTS/Clusters/Methods/w1w124w2626w4309w6550/QM/0DD5496EC209C61A5575F70A3F4669BD
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5068w6550/QM/39BF9404B83932CD3462360F6B7C8A53
$MOLDB_PROJECTS/Clusters/Methods/w1w124w4309w5068w6550/QM/F3B2E3AF8C8DFCEB213D242F2B03A16E
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539/QM/AF55B4A40B81DADB0E13025AB6A0923B
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539/QM/BD47BCDBA8EDBD7599923A46B6AA6D77
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3205w5539w6538/QM/8649CBDF0B3B829CE88D65205195EC89
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3205w5539w6538/QM/621E9ABD14EAF952CAD8BE95BC71166A
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w3205w5539w6538/QM/2B7A45752746069A71D0D279286DE7B5
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w3205w5539w6538/QM/4CFBFFB91E17C16AD101A6484A5ABE13
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w4645w5539w6538/QM/FD4A5898F239BF63F6CF44CEB184E480
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w4645w5539w6538/QM/121520F8EB8DAA431A8231D18FC7E8B1
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539w6538/QM/DD00F376FFD738A25FD2DEA567124392
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5014w5539w6538/QM/3071D8594BE7B733F62D6BB35EE2C45D
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3205w5539w6550/QM/B2CC5005BACD4F005DE6D177A1BDAA8E
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3205w5539w6550/QM/86E08B3E722B61823B8B64186E053188
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3205w5539w6550/QM/A6105B7FDB85E13EFC54D863159DE09D
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3205w5539w6550/QM/97732255E4510C321D6EC0C456191AF6
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5068w5539w6550/QM/A7778936AD9E93BB2BEA2E1A70DA6CF8
$MOLDB_PROJECTS/Clusters/Methods/w1w3205w5068w5539w6550/QM/F3576B7B73FD462893C4387447BD80C5
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w5539w6538/QM/5ECD66701F320744321D174ABB7DB66E
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w5539w6538/QM/8DC40F55A112585714DD508318D33EEF
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w5539w6538/QM/6783F464B4CF32F4725371F361FF3D40
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w5539w6538/QM/17AAB0B7F024CC14DCCCA962454FEA00
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w5539w6538/QM/7E721328C51BDE69A32DF04016C25592
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w5539w6538/QM/5900302B73A90439A99A53F55715113C
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5539w6538/QM/065340B29FA6267F3C0A75F7E7ECD2E5
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5539w6538/QM/079AF367D19A795A86C0E4C7280C8C1A
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w5539w6550/QM/1A2D4D4108D0695195AFC1222D3F46F6
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w5539w6550/QM/DA7020B6FD1DD38D4795EE0CE0C64CBA
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w5539w6550/QM/74E7BD6D2911182D6D2B1C633C1BC472
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w5539w6550/QM/6BA00EBEEAD1E072BF6DDF4F12DA660C
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w5539w6550/QM/3A75949DC60801550EA7559B36619153
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w5539w6550/QM/5342C002B091F1D71CBD18DEB250D94B
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3106w6538/QM/4FD6EE5B318CA205F8391622F5D9D3CD
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w3106w6538/QM/E830CC633F7382EA415B594509723297
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w4645w6538/QM/90DE6EFE8D139BEC5780E8EA983C5DC4
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w4645w6538/QM/6DBA23CD7C622FE9AF2C3CF83C7BD657
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w6538/QM/FA2F54D14CD59F883F1ECF2AA814A609
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5014w6538/QM/486A8D1E469DB3A15EF9A576235555B3
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2293w6538w6550/QM/CC30195FC914F3673EDDC1A80A0BFB42
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2293w6538w6550/QM/0B837CC61BB5BBC7162F6C5A16CC8CB3
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w2626w6538w6550/QM/3E879C6CF60B6A98AB542506B5F4CA7F
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w2626w6538w6550/QM/4A0AFDD41648EB32506C6CD10D340EB8
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5068w6538w6550/QM/661509716AE63E8D03A508B10E866DC1
$MOLDB_PROJECTS/Clusters/Methods/w1w2293w5068w6538w6550/QM/A30BD97560CB0E539527E4D409C950B3
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w4645w6538/QM/F2C74B21A832DD44DF03FF38BDD4C397
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w4645w6538/QM/9A6993A74E2332BE7A6C65CED18603D9
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w6538/QM/A97E8313827192E35161647D7DDB403B
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5014w6538/QM/5329CD3786AD34F4ABFAA639A2215A8C
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3106w6538w6550/QM/D91A22461C4F9CCAB4A5D37F545799A7
$MOLDB_PROJECTS/Clusters/Methods/w1w520w3106w6538w6550/QM/A7BEF39E03D9B4DE96314E3BACC81195
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3106w6538w6550/QM/08CFFC296846DEF0B3F4FF8D56BB625A
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w3106w6538w6550/QM/DE88D5A09222753E43A2F1790C7A26CA
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5068w6538w6550/QM/DB6645D0D91B5433CDED890841CA0643
$MOLDB_PROJECTS/Clusters/Methods/w1w3106w5068w6538w6550/QM/F7DA4D88FE1A2F6DC3086A07C28C7D2A
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w6538/QM/E4B68FF6DA9E78BB9CDC0D78528FF169
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5014w6538/QM/E4F91A8356113FF432D91132B1288259
$MOLDB_PROJECTS/Clusters/Methods/w1w520w4645w6538w6550/QM/CFCF33DE950EAA6CF7BD79DA71251FCF
$MOLDB_PROJECTS/Clusters/Methods/w1w520w4645w6538w6550/QM/EFF5ABE17B6005777A39B87182337C5E
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w4645w6538w6550/QM/E39606740F0C54B8CC88041DA8879EDE
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w4645w6538w6550/QM/28A40C4DDE5343193DA102B41DB4C4A0
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5068w6538w6550/QM/44A511C93957F36372897095E76E1110
$MOLDB_PROJECTS/Clusters/Methods/w1w4645w5068w6538w6550/QM/09B453060B3BCDD079B98CDA8D27E5CE
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w6538w6550/QM/E745EF0B6A88816D28A7A95EA50D19B8
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5014w6538w6550/QM/0BB3953A8EF831DD4E24E16C86D15202
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w6538w6550/QM/462610EC2107181B530319FA8E67DE53
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5014w6538w6550/QM/C5DB6F9839922672BEF23EEED1565AF7
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w6538w6550/QM/0F3C7D5133F2DCDFC3326776ADCEF1EE
$MOLDB_PROJECTS/Clusters/Methods/w1w5014w5068w6538w6550/QM/00DE3583820AA4BBC5455F6B122CACCA
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2626w6550/QM/317871CB9D1505B8E0CFD8F5B5DA9444
$MOLDB_PROJECTS/Clusters/Methods/w1w520w2626w6550/QM/A2F48835B1708461184AAEEE703B27A8
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5068w6550/QM/B5EECB6B41A042DFAA1D5B6F005A0F69
$MOLDB_PROJECTS/Clusters/Methods/w1w520w5068w6550/QM/D1F193D6E24C6A8CD588DD38BF66368B
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5068w6550/QM/725CCEED17D913819E07CC3A895E7C65
$MOLDB_PROJECTS/Clusters/Methods/w1w2626w5068w6550/QM/EB7C32EE5AAED190426A58EDC5F26D46
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
	./prepareInput
	count=$((count+1))
done