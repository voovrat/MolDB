#!/bin/bash
FOLDERS="
$MOLDB_PROJECTS/Clusters/Methods/c1/QM/29A674722D5BAEED665C901C5BC38A60
$MOLDB_PROJECTS/Clusters/Methods/c2/QM/0755597F3C0950EEAA6FDA5683461B03
$MOLDB_PROJECTS/Clusters/Methods/c4/QM/6AE212E53B28B779EE12F038FDAC1A91
$MOLDB_PROJECTS/Clusters/Methods/c5/QM/C52310E1EE77BE27F629109CCCD3493B
$MOLDB_PROJECTS/Clusters/Methods/c0c1/QM/84C219EC72E18387FA5CF37682F008B5
$MOLDB_PROJECTS/Clusters/Methods/c0c2/QM/392E2D79C51233C7EAB3303550C4CC73
$MOLDB_PROJECTS/Clusters/Methods/c0c3/QM/E087B7E58954963EDD021024F072693D
$MOLDB_PROJECTS/Clusters/Methods/c1c3/QM/FB83080011CEA21AAECFACC293CDA98F
$MOLDB_PROJECTS/Clusters/Methods/c0c4/QM/294932F4AD4BDA91691CE5DAB4459BAA
$MOLDB_PROJECTS/Clusters/Methods/c1c4/QM/21C7CD943BD8AA85412A595E341E526D
$MOLDB_PROJECTS/Clusters/Methods/c3c4/QM/B0CE0FF188D4E4146F4C320A836FBD64
$MOLDB_PROJECTS/Clusters/Methods/c0c5/QM/94EEA0727EBF51F256FEC35688E159B3
$MOLDB_PROJECTS/Clusters/Methods/c2c5/QM/71A7A3E25CD4F92D5C6CE83EDA2F2606
$MOLDB_PROJECTS/Clusters/Methods/c3c5/QM/34F41E2EA4F441761161C58F22DBE64E
$MOLDB_PROJECTS/Clusters/Methods/c0c6/QM/DFD48244EBAE79570BC4119463307508
$MOLDB_PROJECTS/Clusters/Methods/c1c6/QM/E06EDA621B3F720330FC5AC362727E4B
$MOLDB_PROJECTS/Clusters/Methods/c3c6/QM/364D591D3301B714C40D3D18D68F3A06
$MOLDB_PROJECTS/Clusters/Methods/c4c6/QM/986C5AC6AF1DA297536985CF5FC6BE13
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2/QM/82CB215467EAAD1F3C3BDB5F8039E53D
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3/QM/9ADAD806F6E0869D57C85A012ABA5121
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3/QM/D1D23BDA832F3AF07225725509083990
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4/QM/CC83E324FA0D38E9818411A6E00CD672
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4/QM/939695361415CAC7BBD37C29BD9E1BFE
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4/QM/EE02D9EF74229FD68F5D18D165AE0046
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4/QM/7F935BFB6E1C8E2AAEE21CD513F46C94
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5/QM/464E0826A5DE9FBE53C2FC45613CFAB4
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5/QM/90289283E8019ABDE55F8232C53604C2
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5/QM/02488074D7B4CBFBF260963D95DC03C9
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5/QM/DF4BDD42759B5BBEFD7A56032EDA0916
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5/QM/F8A9F609B0249761A4A40C8E37131388
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5/QM/B5F510A26BE0F34D2394F962681394E6
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5/QM/B374103262DC90A1A4CF6D0D2A34AF33
$MOLDB_PROJECTS/Clusters/Methods/c0c2c6/QM/E1EA9FBF8F0DB2036CE95C7C55ED73A9
$MOLDB_PROJECTS/Clusters/Methods/c1c2c6/QM/AE753CCC3423AB88EA27F346139A1B9C
$MOLDB_PROJECTS/Clusters/Methods/c1c3c6/QM/2D8265F661DB419F8CAF075676AEB448
$MOLDB_PROJECTS/Clusters/Methods/c2c3c6/QM/7A67F886B8589074F76725E77F421A24
$MOLDB_PROJECTS/Clusters/Methods/c1c4c6/QM/E13F6ED51E0A971529AA5E2136ED27C9
$MOLDB_PROJECTS/Clusters/Methods/c2c4c6/QM/479688194107860D7E5EB4F138CB4337
$MOLDB_PROJECTS/Clusters/Methods/c0c5c6/QM/36342EE37B00E44993AB15F58A578345
$MOLDB_PROJECTS/Clusters/Methods/c1c5c6/QM/458F4B013624E99029AA31DB4B09AFB8
$MOLDB_PROJECTS/Clusters/Methods/c3c5c6/QM/4F52EE90C5C26747155B581C8EBF7EC2
$MOLDB_PROJECTS/Clusters/Methods/c4c5c6/QM/BE25E7C0E57477F6477D97531ED8F8C5
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4/QM/15C010BCB189EEBC745F26620B6B0E6E
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4/QM/C7BFA778CA7A1DB1C9BC14CA902890AA
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4/QM/BF8D011D3FE41CCF18257EA14C41C31B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5/QM/304972031292C5D9885F0E7F2E990ECD
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5/QM/A4D9D8E8DBD1DF7E59FDE70F2F2581C2
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5/QM/47FC8499B9249BCC72E81D199B698ABE
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5/QM/4AE9FFCA79C0C944EBBE127E78A78752
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5/QM/DFA5FC455396EF352E57A5778CA166A8
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5/QM/D4F0CF09ADD03A98E8E24FE41A7E618F
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5/QM/0E08C33A1CD1189BC3E7859C45E1904E
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c6/QM/B7A1380A4723217FAE1146F39C445FC8
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c6/QM/32677F2B668BCAFC8FADBC1562CF983F
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c6/QM/6544550B5A54FD4A3BF990944C6E2C3D
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c6/QM/C0DF2EF031F86D43DFD510D160E09F37
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c6/QM/C3F1AA9120FCE4C2141D6646D243199D
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c6/QM/E1E5C5C1862805BAE018E13F2509EB69
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5c6/QM/CA9B78B8A6F48DDC683981D6E02BBBA1
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5c6/QM/60A15A43655B0499AC923D8A0350D7CD
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5c6/QM/D27340D840BDBECAF4D8308C83DFD2C3
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5c6/QM/F1F806B74D9060AA5160A3E5AEEABB3A
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5c6/QM/B0B2767C4D3B8327810B3BB69D6D123D
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5c6/QM/19088F826B768B275119FC0FF6F7CD00
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5c6/QM/26960CD14AD3AC76255C52EFF22C5755
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4/QM/347E8ACE7073B36E7BF73FDB8F707D0A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5/QM/938AB39FC3A5A038DE9A4E618860154E
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5/QM/2CCAD899BBDE9715582069659C3F1923
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5/QM/0647E228A2813242CAECE6CF9E52195F
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c6/QM/BDC44232D478F8321409E5117ECE9ECE
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c6/QM/66831453635ACF7C9C9DCAD2D6DC19FE
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c6/QM/119C18E677AF7C309A544E86DD96B02A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5c6/QM/BA128281817A167363C88F862387646E
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5c6/QM/DB4FF3E664D816ACDA21485CA0BEDBB2
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5c6/QM/CAD6385E751C99ED0042FCEEB7AD3DAF
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5c6/QM/F39EC2E746F89A23F54AC45D9F017250
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5c6/QM/4F6152AE66CE7CDC1C92C4531DDCA457
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5c6/QM/887D1FD0440FC8D5067C178E8A875B66
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5c6/QM/51401886ECE6261AF1F519757B5C5664
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5/QM/EEC9F83CC75D016FB5ABB32C8D47FF95
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5c6/QM/98CBDA803AB6513598BA87BC68C162AF
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5c6/QM/B6329AE823118D45669E187DE932D2BF
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5c6/QM/C3E3CA78AED5B5D64FF55CA8580AE216
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5c6/QM/466E44593EE5BD5DB752F358F41BDB9A
$MOLDB_PROJECTS/Clusters/Methods/c0/QM/509F7860C3B727CE56FA09AF96F83883
$MOLDB_PROJECTS/Clusters/Methods/c1/QM/C776C8AF7DD326F241A12C2BB169537A
$MOLDB_PROJECTS/Clusters/Methods/c3/QM/B6F4B2A205ECF7309DC563A834420DBD
$MOLDB_PROJECTS/Clusters/Methods/c4/QM/3436BA03B8936E561F155EECC045866B
$MOLDB_PROJECTS/Clusters/Methods/c6/QM/5DE19AB2AD402C6CF9252FEB32DD7794
$MOLDB_PROJECTS/Clusters/Methods/c0c1/QM/9004F5011D6FCE1E471FACAA2F058239
$MOLDB_PROJECTS/Clusters/Methods/c1c2/QM/445F83F9C024011E9EB6EEC2703F0E25
$MOLDB_PROJECTS/Clusters/Methods/c0c3/QM/584E70E523139A10918B3912739CC3CA
$MOLDB_PROJECTS/Clusters/Methods/c2c3/QM/AACA446D9598C8ED5451B9D2E6A9715D
$MOLDB_PROJECTS/Clusters/Methods/c0c4/QM/6E8966D750B5ED728B6ED1796F33E747
$MOLDB_PROJECTS/Clusters/Methods/c2c4/QM/9D50D3961B14ADB6F9505B2586B764FD
$MOLDB_PROJECTS/Clusters/Methods/c3c4/QM/4156D6FD9DA7EBF990755C118A8075DF
$MOLDB_PROJECTS/Clusters/Methods/c1c5/QM/DE2FFDEB0C5B29E7DA5F55A1CC40AB08
$MOLDB_PROJECTS/Clusters/Methods/c2c5/QM/93481BBCEE164D265F1764B017605A9D
$MOLDB_PROJECTS/Clusters/Methods/c4c5/QM/7FE8430BBF2165A622F070ACADF746D0
$MOLDB_PROJECTS/Clusters/Methods/c0c6/QM/8639C36AB0A86DB9590C5093A548504F
$MOLDB_PROJECTS/Clusters/Methods/c2c6/QM/4210FB3B435097C3DE0B05BE997FE1AD
$MOLDB_PROJECTS/Clusters/Methods/c3c6/QM/53768E3427E69F518EF5DE7F940D67F8
$MOLDB_PROJECTS/Clusters/Methods/c5c6/QM/98C1B88F74EC6BD58EAF99A9D87DE324
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2/QM/B2E670D193859CE6207EB2857EB59319
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3/QM/D4773493F0DE29DEC5E742A8662DE5E6
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3/QM/E68E8663EDDF5A5A061325DEAA1DB9AC
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4/QM/2E4DC5A122516371DA69D554BB7C8FF4
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4/QM/431ADAEE0A29E85A0BDD3887E265CA1F
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4/QM/05F4FCF95F1D3A66FF53AFD0B88F515A
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4/QM/98B3558EF91AE1100FC1442CAC40AD51
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5/QM/0CCDF12A2753067CE1AC10F52A146728
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5/QM/7043EF9779CA8BA515D9A7BDC640EF02
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5/QM/2E845DA73155DE71E1416CA1F1F5001D
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5/QM/D935CD1207107AEDA4575B8570A0F080
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5/QM/DC5EFDBBB12D0888D390F808F301CA83
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5/QM/39600F721BEBABDDC3E6581FEE12F1F1
$MOLDB_PROJECTS/Clusters/Methods/c0c1c6/QM/36133A8A5B39D2147572DA9F6BE7B92D
$MOLDB_PROJECTS/Clusters/Methods/c0c2c6/QM/10B82A80E05CAC45666B86CB852930D8
$MOLDB_PROJECTS/Clusters/Methods/c0c3c6/QM/138074C970FA124FAD4D19FE7A76EE93
$MOLDB_PROJECTS/Clusters/Methods/c1c3c6/QM/A2BC7EB0BE351A9DC6E9F5E2EA70FDE1
$MOLDB_PROJECTS/Clusters/Methods/c0c4c6/QM/0116F91753B642D93B37126072D28B76
$MOLDB_PROJECTS/Clusters/Methods/c1c4c6/QM/EAE6B0068144BC4A62A433B087BB7016
$MOLDB_PROJECTS/Clusters/Methods/c3c4c6/QM/27E9B334FD85649051E7E98B2BB1D1CE
$MOLDB_PROJECTS/Clusters/Methods/c0c5c6/QM/B41D1542AEFF28DAA5CB0E0818C893A8
$MOLDB_PROJECTS/Clusters/Methods/c2c5c6/QM/0E59A214EA42BCE6CE4528A274CBAE31
$MOLDB_PROJECTS/Clusters/Methods/c3c5c6/QM/39374924C68B31AE869C01C14A3D1225
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3/QM/D491818CF2B5447A26021923ECA8FE8B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4/QM/340DF1556115384DB6CC2DF29FCBF5DF
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4/QM/00C666C1E369BB4BF38FAEFA94E3B20F
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4/QM/C4109C21CE5FFDF6D4F3BA11CEF60C08
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5/QM/C94F036FE14B172FCD27BA2484581F42
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5/QM/11124ADF5DD120F8521CCBD286C766E8
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5/QM/2B3EF66C0A41CD4B2F86EB2F8FFFB2A9
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5/QM/5AAB90DE2807894AF06C917D9946CB79
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5/QM/C3EC8F5AEA36AE36B63F0B9DF890779D
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5/QM/E0C8F847A642FDDA2DA1EE21DF6A6ADB
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c6/QM/6C0BAC11D2966F0E71198C679486A586
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c6/QM/425B143E2DCF42EE72107500BBAC53A6
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c6/QM/425138BA62719FA11F720DCB8A1017DD
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c6/QM/B3E1466581D3E3856A3C302267C444D2
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c6/QM/F007027BD4FD92FB73AFEC13F93760BE
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c6/QM/EDE533208DC30C994A2078070A2B492D
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c6/QM/F11559BA1C6613A61E81301755832D7F
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5c6/QM/D1B44B3456573E96751D3BE97A461DC0
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5c6/QM/F04F090890F62DFE44006AACA793DF63
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5c6/QM/F874D7824A7D3C20ABA979E7940ABAA5
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5c6/QM/4077B55AD3E8DE1A72E8690A600800BC
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5c6/QM/FB7EB5270396BD90DDFFF7E9B3B69696
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5c6/QM/C33554A276079F3C3065FA5372B1D472
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5c6/QM/36175E92C33FD5B551929B63230897E9
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5/QM/085DE4257181061BC78B824FC6B95FDC
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5/QM/0E33593AB8ADD9A87A971C153673637F
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5/QM/EA0C1279ACBCF9689E3F540166A315AA
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5/QM/202740EBB508365D572DCDAEB1D10AD6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c6/QM/2F810E23ADE508E58BE86931A27C4641
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c6/QM/BAA7F91B6419F7955583DCE339C267D1
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c6/QM/6D0C67990F194523B9DC0CA390F1EAF7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5c6/QM/EFF57746C7E9108E4E1897A8CEC8EADD
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5c6/QM/8B03B902C10670C51F4F50DF8392E208
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5c6/QM/AF0E9013EE457D4058D7DEC288051BA4
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5c6/QM/9734E8E1D9893A45755176A2B7F0E9C0
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5c6/QM/8E6A8AECC542087DA84BFBCFE3D01140
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5c6/QM/A7CE38EFF5062F9B04944420922D1A4C
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5c6/QM/6F38DFAB3FDC03DA163697567825FCDC
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c6/QM/E575355602FD05D917D29A11EC7DB1A7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5c6/QM/0271188CF252EB99FE8EB0B312691649
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5c6/QM/FAF608C72BF3E9FC2B9438EDF7C6925C
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5c6/QM/AD3F4725E0998ADB9315ADB4B21003A6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5c6/QM/9ED6EF7E2C034A98D25ABF87DE84F30E
$MOLDB_PROJECTS/Clusters/Methods/c0/QM/86A0D5C83FA31F0EAA37E3CE8FF0E15A
$MOLDB_PROJECTS/Clusters/Methods/c2/QM/0785CEA35CA247229D54A91C7101E1EC
$MOLDB_PROJECTS/Clusters/Methods/c3/QM/1164C59A5818F6CBFE92F34F540380BE
$MOLDB_PROJECTS/Clusters/Methods/c5/QM/749CF9FCAAA23791355A26D9EA9C5E71
$MOLDB_PROJECTS/Clusters/Methods/c6/QM/DAE6E42911362434F81CBC16E5C683EB
$MOLDB_PROJECTS/Clusters/Methods/c0c2/QM/BA3626DF7CC44D7E070E914A32FA3E92
$MOLDB_PROJECTS/Clusters/Methods/c1c2/QM/CEBDD1DD6F3637DB6353EC9323265474
$MOLDB_PROJECTS/Clusters/Methods/c1c3/QM/0B9789CE1AD873DED8C9DA0D6C04EE1E
$MOLDB_PROJECTS/Clusters/Methods/c2c3/QM/A6BCE09B51465BCCE8F083515D83BBF0
$MOLDB_PROJECTS/Clusters/Methods/c1c4/QM/177FD008C8A492FA60E4790A5BB50BC6
$MOLDB_PROJECTS/Clusters/Methods/c2c4/QM/2B4B697494A1DC76D6A661E778A79F45
$MOLDB_PROJECTS/Clusters/Methods/c0c5/QM/756D00F593B589430FD6EEA3501BF6EB
$MOLDB_PROJECTS/Clusters/Methods/c1c5/QM/8E8D404EB503795FF1C98BAB302120F5
$MOLDB_PROJECTS/Clusters/Methods/c3c5/QM/ADB93DA249C649AE7952C9308E5A1C2C
$MOLDB_PROJECTS/Clusters/Methods/c4c5/QM/D58554E5150F2AD2E44DA390435E98E3
$MOLDB_PROJECTS/Clusters/Methods/c1c6/QM/7325E3B82E3E631D0D636725AB819C8D
$MOLDB_PROJECTS/Clusters/Methods/c2c6/QM/00EEB84CD260DDF0676EE71CC5F7C2D1
$MOLDB_PROJECTS/Clusters/Methods/c4c6/QM/4073A60272DEFD4817EFB3241A1F4A44
$MOLDB_PROJECTS/Clusters/Methods/c5c6/QM/EC64365D9CE9C6EC7CE3C2754847717C
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3/QM/552989EB3626A51BF19847B162008B29
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3/QM/E1E7BED4D593799A1968CFE30E2B9245
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4/QM/A5883F8AB48B88D927E37C5AB1584339
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4/QM/5BD08CD7086CBECC04BC7C80FE053DAC
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4/QM/1990303903953FCF16F08E58BB0E3018
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4/QM/8649B38CA250B9EFA11B9AFB4033CF64
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5/QM/E0EA692CE9A8A5EE6ACC9A368E531D7E
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5/QM/C1320779FB3D93517F174A139DAFA6FE
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5/QM/05E081A992AD29486A803AAC37680100
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5/QM/5EF17D3E52CEB915F760AD47502AA81B
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5/QM/D7978A4F2167806FB2EAE597F6138218
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5/QM/FC99E7786D3012541385A45436A889A0
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5/QM/B19F9F250E5C6A6BA1901917E4B758A3
$MOLDB_PROJECTS/Clusters/Methods/c0c1c6/QM/FE2B562B37FC5AD2EDCAB356B9B585B6
$MOLDB_PROJECTS/Clusters/Methods/c1c2c6/QM/14FE04F611A168D63ECAEC80908F9A30
$MOLDB_PROJECTS/Clusters/Methods/c0c3c6/QM/00691DC6F49B72BDFABD8ADB0BDF552C
$MOLDB_PROJECTS/Clusters/Methods/c2c3c6/QM/8FEB118C4FCF971ABE0573D3359C8B60
$MOLDB_PROJECTS/Clusters/Methods/c0c4c6/QM/D3AE3DECA218EC950EC2BC1CD72A43D2
$MOLDB_PROJECTS/Clusters/Methods/c2c4c6/QM/A1DBA8775CA2C89CEED78E4E37BA4F01
$MOLDB_PROJECTS/Clusters/Methods/c3c4c6/QM/94356AA1F132AA0B568A872CEC1F7C07
$MOLDB_PROJECTS/Clusters/Methods/c1c5c6/QM/770EA90B793909EBD625EAE914ED98B1
$MOLDB_PROJECTS/Clusters/Methods/c2c5c6/QM/E7DD6F55DB28A49B3325E68F2F716319
$MOLDB_PROJECTS/Clusters/Methods/c4c5c6/QM/899A458D1ED51DD32EDBA6149C21BB6A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3/QM/6BB724B50F1BE9EE06A93DC969B460C6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4/QM/54339F7AE0C2BCD7FFAFF0C490AFB3B2
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4/QM/AF2BCEE14D4F9A587F22AC9788C1A756
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5/QM/EB41BD464412C4F81B82474E1F52BC5B
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5/QM/B9EBBF7A3C4A1F0905AA3B208FC6F9D5
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5/QM/4CA0D602926478AC6F1C12B651A427A6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5/QM/C80CDAAAABB332F8F5A82A9AF93C61FD
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5/QM/836ADB6087A920E023CD7C1DB7C2EC6A
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5/QM/A84A1B5DFB8D3FA785A1252E16F02436
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5/QM/651F6C25B66EDC8722FA586AA17ACF17
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c6/QM/99BA56EF96C2EAD99805EC1A347190FB
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c6/QM/5A5BB93B63C885D550A6C56FA51333D2
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c6/QM/A164D6CC3A81DDF59D4E3BBE505AA35B
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c6/QM/0B674EAC56336091858636D4E9185DEE
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c6/QM/9765F54FA1B75F25FA5BCD3A827EBDB6
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c6/QM/E447F88FC7577A6FF928113CA9ED9257
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c6/QM/88B32CC50E9534962AD1564F13B9B8BE
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5c6/QM/C093618793A3915EF041667FF07C2403
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5c6/QM/43E9EF3CFB6F881C179355C42C259195
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5c6/QM/1D31FF25126721B835A6F5C03160BEF5
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5c6/QM/030DDDCCC945350F476F4C3C46C5E8E3
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5c6/QM/B4FEE05E5B92E7566F46044FC8482E46
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5c6/QM/FC82E5F791C2EAD5A7DC012DD6E569A3
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4/QM/5A074DEB4C3FE61B673EF80223C9BFF6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5/QM/2BDD9121D03275E8D9C74F6A07E9A37A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5/QM/482614574A4D44AEA605857B76881EC3
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5/QM/BF7CF32A36B67982D6970D5E1414A487
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c6/QM/62776EE29243620DA348CD2D04B7108D
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c6/QM/53CA017894A225B2D3858D8444D3BE49
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c6/QM/B1E073424D4D1EE33A830ACE332070DC
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c6/QM/3353229038FBC5D36ED8E424522D86E8
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5c6/QM/5E3C8F6F309F3C6D83DFF57B93B5BF45
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5c6/QM/23234FB1186956FA5431DC4065546F02
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5c6/QM/19F5EC51DFCE2221A9889C00BDDDF418
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5c6/QM/51644C914A51D926EFF92E8AC4B7E9C6
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5c6/QM/685FE46F60D7A717823B07118A21C7BD
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5c6/QM/287584444031ECC2318B8144C3E5BF38
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5/QM/E3DF0F86D154A88916CE91C37BD2B918
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c6/QM/F6A47724AE367497D1A1926F15712DB7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5c6/QM/7CBD4B79E015B5D6958745492B45E735
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5c6/QM/CD06EC23596490FC8268EC822EF4E5CA
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5c6/QM/AC3085613C047D5BC32AD32565559F01
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5c6/QM/C17BB5AE9AA2F479376F6165849356C5
$MOLDB_PROJECTS/Clusters/Methods/c1/QM/DBBABF7307B6AF253CC32D3760A36027
$MOLDB_PROJECTS/Clusters/Methods/c2/QM/BF8E41406E5E015B3238CE660F0E3AD7
$MOLDB_PROJECTS/Clusters/Methods/c4/QM/7B956D2D677814B5E7DE35ADF9B96DEE
$MOLDB_PROJECTS/Clusters/Methods/c5/QM/25C8ED8E2F188D0D562F25205E7AFAE3
$MOLDB_PROJECTS/Clusters/Methods/c0c1/QM/BE146582ECC7CA1E740AB0098BD9CF6C
$MOLDB_PROJECTS/Clusters/Methods/c0c2/QM/50AF57FCF0DDB60006BB4E355CEB1BF7
$MOLDB_PROJECTS/Clusters/Methods/c0c3/QM/DA77930E1E96154DC1F84D5750FCDC1C
$MOLDB_PROJECTS/Clusters/Methods/c1c3/QM/0482F66C5E8A69A5DCE3F7163CBFCAD6
$MOLDB_PROJECTS/Clusters/Methods/c0c4/QM/525B208947C3CE6829E7061BE88D38E8
$MOLDB_PROJECTS/Clusters/Methods/c1c4/QM/23D1E5AA41E4BBFD827A4754EBF3367E
$MOLDB_PROJECTS/Clusters/Methods/c3c4/QM/AFE19E5DBA0630A72665D3F70E24C77F
$MOLDB_PROJECTS/Clusters/Methods/c0c5/QM/1F6683B46D9C81F90370A9DBA4936C20
$MOLDB_PROJECTS/Clusters/Methods/c2c5/QM/3E0B1A927A1766AF802BBD8AFFE015F0
$MOLDB_PROJECTS/Clusters/Methods/c3c5/QM/A99ABEDF97579B8554B999EEF47EBE09
$MOLDB_PROJECTS/Clusters/Methods/c0c6/QM/AA9D0446331F5845FBFF4C05EB0CAF0A
$MOLDB_PROJECTS/Clusters/Methods/c1c6/QM/2496A3F7173EC17B1C8A264CB108810C
$MOLDB_PROJECTS/Clusters/Methods/c3c6/QM/E3C5043D8651B6CA3485C30B2FC08FE7
$MOLDB_PROJECTS/Clusters/Methods/c4c6/QM/5177F83E084432C6CD30874EEE68545A
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2/QM/043BD6050AA85C5B44DE905A5881D0E8
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3/QM/99E48D49B8C451522C00D2F51A31D7E4
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3/QM/3210EC3A8909F9939F2300C8854F18EE
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4/QM/3292F4CF47E85A0A925DDF34A64AB411
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4/QM/0F5B713EFC847028F7EEC21FC5D5FA35
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4/QM/16C53A70A1A23148EF91298D733E86AA
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4/QM/C47EB1A5A8DC8F384588AD5DBF4353E7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5/QM/DC1A0CA6E3C5013E2D66165CD59D877E
$MOLDB_PROJECTS/Clusters/Methods/c1c2c5/QM/228E609050636E6ADBB1038FB55BB4C6
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5/QM/FA0B0307C5D32F2915A08324E9F8E893
$MOLDB_PROJECTS/Clusters/Methods/c2c3c5/QM/17D656F45E5D37F03A7D7707E1C8C405
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5/QM/5520CD9FAECCE1D578976612494A715D
$MOLDB_PROJECTS/Clusters/Methods/c2c4c5/QM/A0C032769F54E7928CEE6DD47341DF30
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5/QM/3A75B207A7BD977C988A0CDF1CCAD9FF
$MOLDB_PROJECTS/Clusters/Methods/c0c2c6/QM/97D845A41A36A3A1883F42322606761C
$MOLDB_PROJECTS/Clusters/Methods/c1c2c6/QM/E93FF4AF0314ADEB5AF5ACB847C398C1
$MOLDB_PROJECTS/Clusters/Methods/c1c3c6/QM/9DC4C6C297F3CD7FEFB248B0D17563D1
$MOLDB_PROJECTS/Clusters/Methods/c2c3c6/QM/BEE69372004BBB6BD967195B179B7C3A
$MOLDB_PROJECTS/Clusters/Methods/c1c4c6/QM/9C9FB3ECA9FFDBF316AF58F19CBB5983
$MOLDB_PROJECTS/Clusters/Methods/c2c4c6/QM/5371EC3E6636A951912E885B89B435C7
$MOLDB_PROJECTS/Clusters/Methods/c0c5c6/QM/43C38ABC044635D7A7BDF084E6EF09FC
$MOLDB_PROJECTS/Clusters/Methods/c1c5c6/QM/F12528D7B7A6AC14F70E516CA6EB2E31
$MOLDB_PROJECTS/Clusters/Methods/c3c5c6/QM/6BB64E59931B60EF4E933EBE6F3C54BB
$MOLDB_PROJECTS/Clusters/Methods/c4c5c6/QM/299FB0D6A1A2E852DF165097F3CC3A63
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4/QM/517C5755D300C5370D7FB02380190017
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4/QM/E8431E7DA6F9D4996BB5CA95DD118EE0
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4/QM/71E07B9C9CDE6ED0E20F6E24035F3B96
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5/QM/AF8AD112CA27307807560A30EFDA4640
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c5/QM/4E39182EE84B255666C348510C4FECB7
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5/QM/2FEA4F02BAC10533D2034633C13D814A
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c5/QM/D5E64BBB37C933489882C24B340DCED7
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5/QM/E9DD76434F87988DF7193AA269113A29
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c5/QM/E4249439389122FECA4536BA2A440CE3
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5/QM/3C230B5662936DF491665988B01E7FC8
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c6/QM/C96E5CD9533226BD2B0AC047865BEDDE
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c6/QM/2158A7034AF8743A368CA0C310E72980
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c6/QM/FDCF4B856055E0A1C37C72149D7035BD
$MOLDB_PROJECTS/Clusters/Methods/c0c2c4c6/QM/217B5ACDE62E1B54EB4D9C44658A1385
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c6/QM/8277F2BA073347837D7E253D3153F811
$MOLDB_PROJECTS/Clusters/Methods/c1c3c4c6/QM/87D196C8AE82F7E420AD5690EDBB7911
$MOLDB_PROJECTS/Clusters/Methods/c0c1c5c6/QM/1F59478AA939ECEDD4902AF977A26618
$MOLDB_PROJECTS/Clusters/Methods/c0c2c5c6/QM/93B59E71980155CB3334B6563DD4B3DB
$MOLDB_PROJECTS/Clusters/Methods/c0c3c5c6/QM/E62531CED75D1C0CADD847A805516D99
$MOLDB_PROJECTS/Clusters/Methods/c1c3c5c6/QM/7EB8BA2824AF3F2921532FF9E0BCB8AF
$MOLDB_PROJECTS/Clusters/Methods/c0c4c5c6/QM/841AC9742D9E0AD51B2463337AEB3935
$MOLDB_PROJECTS/Clusters/Methods/c1c4c5c6/QM/5C94135B7922545DA9D799C51EA0C797
$MOLDB_PROJECTS/Clusters/Methods/c3c4c5c6/QM/9960F20666F3A2C76B722AC302188FD7
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4/QM/1125CE77311C478A3400FE996BD671B6
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5/QM/90F9B2F7B56FE4F82F62AF719F24A5BA
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c5/QM/FBC2418FFCE8147BA4531CE79E459CAB
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5/QM/F4574055030A5637069B5C39DD4C50F4
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c6/QM/3DAD7CC7655179EAFADF08E0C514E2FC
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c4c6/QM/EF870993023CAA21B08311067BC8BAE1
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c6/QM/F6F32F25617EDA1396C1A6440D906978
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c5c6/QM/94FD14539871C0868C0077915F10BF21
$MOLDB_PROJECTS/Clusters/Methods/c0c1c3c5c6/QM/47A9925FD6FFC05BC27A3C820B1CB0DC
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c5c6/QM/0B9F3FE255062E64843B2419EA4916A9
$MOLDB_PROJECTS/Clusters/Methods/c0c1c4c5c6/QM/5FC4E04871AC311F016C9F4522D5621D
$MOLDB_PROJECTS/Clusters/Methods/c1c2c4c5c6/QM/A9A040E29A1B9807000D2E8A8300C97C
$MOLDB_PROJECTS/Clusters/Methods/c0c3c4c5c6/QM/7473C36AC0ED697DAE97F719015D59B5
$MOLDB_PROJECTS/Clusters/Methods/c2c3c4c5c6/QM/ABDBF2CDB862EA3DECC83A1DC70AB157
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c4c5/QM/EC2823463CAC4DB6B7E85390833A7300
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c3c5c6/QM/D4EEBD1052F2738ADE0CB94E4C0B1411
$MOLDB_PROJECTS/Clusters/Methods/c0c1c2c4c5c6/QM/F2F29C1A69DED6DE164B36AA33662693
$MOLDB_PROJECTS/Clusters/Methods/c0c2c3c4c5c6/QM/8F6F7AFBD06B66323B099E8FAC38757D
$MOLDB_PROJECTS/Clusters/Methods/c1c2c3c4c5c6/QM/F6ADF659670A17E970E0D691E5154A1E
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
