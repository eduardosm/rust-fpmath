use crate::double::Double;
use crate::scalbn;
use crate::traits::Float as _;

/// Returns `(k, lo, ln_hi)` such as `ln(x) = k * ln(2) + ln_hi + ln(lo)`
#[inline]
pub(crate) fn log_core(x: Double<f64>) -> (f64, Double<f64>, Double<f64>) {
    // GENERATE: ln_table Double<f64> 6
    // LN_TBL[i] = ln(1 + i / 64)       if i < 32
    //           = ln((1 + i / 64) / 2) if i >= 32
    static LN_TBL: [Double<f64>; 64] = [
        Double::new(
            f64::from_bits(0x0000000000000000),
            f64::from_bits(0x0000000000000000),
        ), // 0
        Double::new(
            f64::from_bits(0x3F8FC0A8B0FC03E3),
            f64::from_bits(0x3C39F3DB4E9A6F58),
        ), // 1.550418653596525415085404604245e-2
        Double::new(
            f64::from_bits(0x3F9F829B0E783300),
            f64::from_bits(0x3C333E3F04F1EF23),
        ), // 3.077165866675368837102820759677e-2
        Double::new(
            f64::from_bits(0x3FA77458F632DCFC),
            f64::from_bits(0x3C418D3CA87B9296),
        ), // 4.580953603129420316667926761466e-2
        Double::new(
            f64::from_bits(0x3FAF0A30C01162A6),
            f64::from_bits(0x3C485F325C5BBACD),
        ), // 6.062462181643484258060613204042e-2
        Double::new(
            f64::from_bits(0x3FB341D7961BD1D0),
            f64::from_bits(0x3C6253306EC209A2),
        ), // 7.522342123758752569860533998366e-2
        Double::new(
            f64::from_bits(0x3FB6F0D28AE56B4B),
            f64::from_bits(0x3C637C93373DA337),
        ), // 8.961215868968713261995146937848e-2
        Double::new(
            f64::from_bits(0x3FBA926D3A4AD563),
            f64::from_bits(0x3C5942F48AA70EA9),
        ), // 1.037967936816435648260618037640e-1
        Double::new(
            f64::from_bits(0x3FBE27076E2AF2E5),
            f64::from_bits(0x3C6D3D50FFFC3FD4),
        ), // 1.177830356563834545387941094705e-1
        Double::new(
            f64::from_bits(0x3FC0D77E7CD08E59),
            f64::from_bits(0x3C69A5DC5E9030AC),
        ), // 1.315763577887192725887161286895e-1
        Double::new(
            f64::from_bits(0x3FC29552F81FF523),
            f64::from_bits(0x3C6301771C407DBF),
        ), // 1.451820098444978972819350637406e-1
        Double::new(
            f64::from_bits(0x3FC44D2B6CCB7D1E),
            f64::from_bits(0x3C69F4F6543E1F88),
        ), // 1.586050301766385840933711746258e-1
        Double::new(
            f64::from_bits(0x3FC5FF3070A793D3),
            f64::from_bits(0x3C790E7C4140E424),
        ), // 1.718502569266592223400989460551e-1
        Double::new(
            f64::from_bits(0x3FC7AB890210D909),
            f64::from_bits(0x3C4BE36B2D6A0608),
        ), // 1.849223384940119926639035926592e-1
        Double::new(
            f64::from_bits(0x3FC9525A9CF456B4),
            f64::from_bits(0x3C6D904C1D4E2E26),
        ), // 1.978257433299198803625720711970e-1
        Double::new(
            f64::from_bits(0x3FCAF3C94E80BFF2),
            f64::from_bits(0x3C7B19CC0326F99F),
        ), // 2.105647691073496376695528127324e-1
        Double::new(
            f64::from_bits(0x3FCC8FF7C79A9A21),
            f64::from_bits(0x3C7584BB03DE5FF7),
        ), // 2.231435513142097557662950903098e-1
        Double::new(
            f64::from_bits(0x3FCE27076E2AF2E5),
            f64::from_bits(0x3C7D3D50FFFC3FD4),
        ), // 2.355660713127669090775882189410e-1
        Double::new(
            f64::from_bits(0x3FCFB9186D5E3E2A),
            f64::from_bits(0x3C71AAA8CD86F29A),
        ), // 2.478361639045812567806027657465e-1
        Double::new(
            f64::from_bits(0x3FD0A324E27390E3),
            f64::from_bits(0x3C77DCFDE8061C03),
        ), // 2.599575244369260669720794945423e-1
        Double::new(
            f64::from_bits(0x3FD1675CABABA60E),
            f64::from_bits(0x3C2CE63EAB883717),
        ), // 2.719337154836417588316694945330e-1
        Double::new(
            f64::from_bits(0x3FD22941FBCF7965),
            f64::from_bits(0x3C844850A7B4EBA8),
        ), // 2.837681731306445983469012223503e-1
        Double::new(
            f64::from_bits(0x3FD2E8E2BAE11D30),
            f64::from_bits(0x3C8385992350A103),
        ), // 2.954642128938358763866819060550e-1
        Double::new(
            f64::from_bits(0x3FD3A64C556945E9),
            f64::from_bits(0x3C88E5E6B9AE81A1),
        ), // 3.070250352949118620751245405354e-1
        Double::new(
            f64::from_bits(0x3FD4618BC21C5EC2),
            f64::from_bits(0x3C7F42DECDECCF1D),
        ), // 3.184537311185346158102472135906e-1
        Double::new(
            f64::from_bits(0x3FD51AAD872DF82D),
            f64::from_bits(0x3C43927AC19F55E3),
        ), // 3.297532863724679818144228119208e-1
        Double::new(
            f64::from_bits(0x3FD5D1BDBF5809CA),
            f64::from_bits(0x3C74236383DC7FE1),
        ), // 3.409265869705932103050891997804e-1
        Double::new(
            f64::from_bits(0x3FD686C81E9B14AE),
            f64::from_bits(0x3C888857C2029C71),
        ), // 3.519764231571781846554474562594e-1
        Double::new(
            f64::from_bits(0x3FD739D7F6BBD006),
            f64::from_bits(0x3C839C498A7F5A7E),
        ), // 3.629054936893684531378243459775e-1
        Double::new(
            f64::from_bits(0x3FD7EAF83B82AFC3),
            f64::from_bits(0x3C792CE979ED2950),
        ), // 3.737164097935840808210168327158e-1
        Double::new(
            f64::from_bits(0x3FD89A3386C1425A),
            f64::from_bits(0x3C86B4E310220783),
        ), // 3.844116989103320397347900624813e-1
        Double::new(
            f64::from_bits(0x3FD947941C2116FA),
            f64::from_bits(0x3C87499BA28FA20E),
        ), // 3.949938082408689781063940363650e-1
        Double::new(
            f64::from_bits(0xBFD269621134DB92),
            f64::from_bits(0xBC7E0EFADD9DB02B),
        ), // -2.876820724517809274392190059938e-1
        Double::new(
            f64::from_bits(0xBFD1BF99635A6B94),
            f64::from_bits(0xBC8BB5451EF6DB77),
        ), // -2.773192854162343438039032285033e-1
        Double::new(
            f64::from_bits(0xBFD1178E8227E47B),
            f64::from_bits(0xBC8BC671683F8E5C),
        ), // -2.670627852490452462926872418627e-1
        Double::new(
            f64::from_bits(0xBFD07138604D5862),
            f64::from_bits(0xBC7CDB16ED4E9138),
        ), // -2.569104137850272390681907983971e-1
        Double::new(
            f64::from_bits(0xBFCF991C6CB3B379),
            f64::from_bits(0xBC6F665066F980A2),
        ), // -2.468600779315257978846419408385e-1
        Double::new(
            f64::from_bits(0xBFCE530EFFE71012),
            f64::from_bits(0xBC42276041F43042),
        ), // -2.369097470783577150364265832942e-1
        Double::new(
            f64::from_bits(0xBFCD1037F2655E7B),
            f64::from_bits(0xBC660629242471A2),
        ), // -2.270574506353460848586128739534e-1
        Double::new(
            f64::from_bits(0xBFCBD087383BD8AD),
            f64::from_bits(0xBC3DD355F6A516D7),
        ), // -2.173012756899813951520225351538e-1
        Double::new(
            f64::from_bits(0xBFCA93ED3C8AD9E3),
            f64::from_bits(0xBC6BCAFA9DE97203),
        ), // -2.076393647782445016154410442674e-1
        Double::new(
            f64::from_bits(0xBFC95A5ADCF7017F),
            f64::from_bits(0xBC5142C507FB7A3D),
        ), // -1.980699137620937948192675366153e-1
        Double::new(
            f64::from_bits(0xBFC823C16551A3C1),
            f64::from_bits(0xBC776E698C7A0C3F),
        ), // -1.885911698075500223589235897200e-1
        Double::new(
            f64::from_bits(0xBFC6F0128B756ABB),
            f64::from_bits(0xBC7390D31EF0F4CA),
        ), // -1.792014294577109926162260331840e-1
        Double::new(
            f64::from_bits(0xBFC5BF406B543DB1),
            f64::from_bits(0xBC7F70525D9F9041),
        ), // -1.698990367953974729004248965233e-1
        Double::new(
            f64::from_bits(0xBFC4913D8333B560),
            f64::from_bits(0xBC7BCAA7EDB3C3B3),
        ), // -1.606823816904734655433083979980e-1
        Double::new(
            f64::from_bits(0xBFC365FCB0159016),
            f64::from_bits(0xBC57D411A5B944AD),
        ), // -1.515498981272009378406898175577e-1
        Double::new(
            f64::from_bits(0xBFC23D712A49C201),
            f64::from_bits(0xBC748E3F4F7D714B),
        ), // -1.425000626072830301572839422533e-1
        Double::new(
            f64::from_bits(0xBFC1178E8227E47B),
            f64::from_bits(0xBC7BC671683F8E5C),
        ), // -1.335313926245226231463436209313e-1
        Double::new(
            f64::from_bits(0xBFBFE89139DBD565),
            f64::from_bits(0xBC629B05EF503636),
        ), // -1.246424452072765973384933565912e-1
        Double::new(
            f64::from_bits(0xBFBDA727638446A2),
            f64::from_bits(0xBC5401FA71733019),
        ), // -1.158318155251217050991200599387e-1
        Double::new(
            f64::from_bits(0xBFBB6AC88DAD5B1B),
            f64::from_bits(0xBC6BFEA044B8D698),
        ), // -1.070981355563671005131126851709e-1
        Double::new(
            f64::from_bits(0xBFB9335E5D594988),
            f64::from_bits(0xBC65C3ABD47D99A5),
        ), // -9.844007281325251990288857492897e-2
        Double::new(
            f64::from_bits(0xBFB700D30AEAC0E0),
            f64::from_bits(0xBC6E8DA99DED3230),
        ), // -8.985632912186104707664693479687e-2
        Double::new(
            f64::from_bits(0xBFB4D3115D207EAC),
            f64::from_bits(0xBC5769F42C7842CC),
        ), // -8.134563945395240588734235502936e-2
        Double::new(
            f64::from_bits(0xBFB2AA04A44717A4),
            f64::from_bits(0xBC61751639682E04),
        ), // -7.290677080808778056573748889093e-2
        Double::new(
            f64::from_bits(0xBFB08598B59E3A06),
            f64::from_bits(0xBC61147FB37EA067),
        ), // -6.453852113757117167292391568399e-2
        Double::new(
            f64::from_bits(0xBFACCB73CDDDB2CB),
            f64::from_bits(0xBC50DB827D7F8816),
        ), // -5.623971832287607779673769427698e-2
        Double::new(
            f64::from_bits(0xBFA894AA149FB343),
            f64::from_bits(0xBC3A8BE97660A23D),
        ), // -4.800921918636060775200362532344e-2
        Double::new(
            f64::from_bits(0xBFA466AED42DE3E9),
            f64::from_bits(0xBC519148405AF641),
        ), // -3.984590854719967065861624024730e-2
        Double::new(
            f64::from_bits(0xBFA0415D89E74444),
            f64::from_bits(0xBC4C05CF1D753622),
        ), // -3.174869831458030115699628274853e-2
        Double::new(
            f64::from_bits(0xBF98492528C8CABE),
            f64::from_bits(0xBC4173697CF302CD),
        ), // -2.371652661731604211834685052867e-2
        Double::new(
            f64::from_bits(0xBF90205658935847),
            f64::from_bits(0xBC327C8E8416E71F),
        ), // -1.574835696813916860754951146083e-2
        Double::new(
            f64::from_bits(0xBF8010157588DE71),
            f64::from_bits(0xBC146662D417CED0),
        ), // -7.843177461025892873184042490944e-3
    ];

    // GENERATE: ln_lo_scale_table Double<f64> 6
    // LN_LO_SCALE_TBL[i] = 1 / (1 + i / 64)
    static LN_LO_SCALE_TBL: [Double<f64>; 64] = [
        Double::new(
            f64::from_bits(0x3FF0000000000000),
            f64::from_bits(0x0000000000000000),
        ), // 1.000000000000000000000000000000e0
        Double::new(
            f64::from_bits(0x3FEF81F81F81F81F),
            f64::from_bits(0x3C903F03F03F03F0),
        ), // 9.846153846153846153846153846154e-1
        Double::new(
            f64::from_bits(0x3FEF07C1F07C1F07),
            f64::from_bits(0x3C983E0F83E0F83E),
        ), // 9.696969696969696969696969696970e-1
        Double::new(
            f64::from_bits(0x3FEE9131ABF0B767),
            f64::from_bits(0x3C7503D226357E17),
        ), // 9.552238805970149253731343283582e-1
        Double::new(
            f64::from_bits(0x3FEE1E1E1E1E1E1E),
            f64::from_bits(0x3C6E1E1E1E1E1E1E),
        ), // 9.411764705882352941176470588235e-1
        Double::new(
            f64::from_bits(0x3FEDAE6076B981DA),
            f64::from_bits(0x3C9CC0ED7303B5CC),
        ), // 9.275362318840579710144927536232e-1
        Double::new(
            f64::from_bits(0x3FED41D41D41D41D),
            f64::from_bits(0x3C80750750750750),
        ), // 9.142857142857142857142857142857e-1
        Double::new(
            f64::from_bits(0x3FECD85689039B0A),
            f64::from_bits(0x3C9A240E6C2B4482),
        ), // 9.014084507042253521126760563380e-1
        Double::new(
            f64::from_bits(0x3FEC71C71C71C71C),
            f64::from_bits(0x3C8C71C71C71C71C),
        ), // 8.888888888888888888888888888889e-1
        Double::new(
            f64::from_bits(0x3FEC0E070381C0E0),
            f64::from_bits(0x3C8C0E070381C0E0),
        ), // 8.767123287671232876712328767123e-1
        Double::new(
            f64::from_bits(0x3FEBACF914C1BACF),
            f64::from_bits(0x3C922983759F2298),
        ), // 8.648648648648648648648648648649e-1
        Double::new(
            f64::from_bits(0x3FEB4E81B4E81B4E),
            f64::from_bits(0x3C90369D0369D037),
        ), // 8.533333333333333333333333333333e-1
        Double::new(
            f64::from_bits(0x3FEAF286BCA1AF28),
            f64::from_bits(0x3C8AF286BCA1AF28),
        ), // 8.421052631578947368421052631579e-1
        Double::new(
            f64::from_bits(0x3FEA98EF606A63BD),
            f64::from_bits(0x3C903531DEC0D4C7),
        ), // 8.311688311688311688311688311688e-1
        Double::new(
            f64::from_bits(0x3FEA41A41A41A41A),
            f64::from_bits(0x3C80690690690690),
        ), // 8.205128205128205128205128205128e-1
        Double::new(
            f64::from_bits(0x3FE9EC8E951033D9),
            f64::from_bits(0x3C6D2A2067B23A54),
        ), // 8.101265822784810126582278481013e-1
        Double::new(
            f64::from_bits(0x3FE9999999999999),
            f64::from_bits(0x3C93333333333333),
        ), // 8.000000000000000000000000000000e-1
        Double::new(
            f64::from_bits(0x3FE948B0FCD6E9E0),
            f64::from_bits(0x3C8948B0FCD6E9E0),
        ), // 7.901234567901234567901234567901e-1
        Double::new(
            f64::from_bits(0x3FE8F9C18F9C18F9),
            f64::from_bits(0x3C9831F3831F3832),
        ), // 7.804878048780487804878048780488e-1
        Double::new(
            f64::from_bits(0x3FE8ACB90F6BF3A9),
            f64::from_bits(0x3C946F0940C565C8),
        ), // 7.710843373493975903614457831325e-1
        Double::new(
            f64::from_bits(0x3FE8618618618618),
            f64::from_bits(0x3C88618618618618),
        ), // 7.619047619047619047619047619048e-1
        Double::new(
            f64::from_bits(0x3FE8181818181818),
            f64::from_bits(0x3C68181818181818),
        ), // 7.529411764705882352941176470588e-1
        Double::new(
            f64::from_bits(0x3FE7D05F417D05F4),
            f64::from_bits(0x3C67D05F417D05F4),
        ), // 7.441860465116279069767441860465e-1
        Double::new(
            f64::from_bits(0x3FE78A4C8178A4C8),
            f64::from_bits(0x3C678A4C8178A4C8),
        ), // 7.356321839080459770114942528736e-1
        Double::new(
            f64::from_bits(0x3FE745D1745D1745),
            f64::from_bits(0x3C9A2E8BA2E8BA2F),
        ), // 7.272727272727272727272727272727e-1
        Double::new(
            f64::from_bits(0x3FE702E05C0B8170),
            f64::from_bits(0x3C7702E05C0B8170),
        ), // 7.191011235955056179775280898876e-1
        Double::new(
            f64::from_bits(0x3FE6C16C16C16C16),
            f64::from_bits(0x3C982D82D82D82D8),
        ), // 7.111111111111111111111111111111e-1
        Double::new(
            f64::from_bits(0x3FE6816816816816),
            f64::from_bits(0x3C902D02D02D02D0),
        ), // 7.032967032967032967032967032967e-1
        Double::new(
            f64::from_bits(0x3FE642C8590B2164),
            f64::from_bits(0x3C7642C8590B2164),
        ), // 6.956521739130434782608695652174e-1
        Double::new(
            f64::from_bits(0x3FE6058160581605),
            f64::from_bits(0x3C902C0B02C0B02C),
        ), // 6.881720430107526881720430107527e-1
        Double::new(
            f64::from_bits(0x3FE5C9882B931057),
            f64::from_bits(0x3C7310572620AE4C),
        ), // 6.808510638297872340425531914894e-1
        Double::new(
            f64::from_bits(0x3FE58ED2308158ED),
            f64::from_bits(0x3C71840AC7691841),
        ), // 6.736842105263157894736842105263e-1
        Double::new(
            f64::from_bits(0x3FE5555555555555),
            f64::from_bits(0x3C85555555555555),
        ), // 6.666666666666666666666666666667e-1
        Double::new(
            f64::from_bits(0x3FE51D07EAE2F815),
            f64::from_bits(0x3C6D07EAE2F8151D),
        ), // 6.597938144329896907216494845361e-1
        Double::new(
            f64::from_bits(0x3FE4E5E0A72F0539),
            f64::from_bits(0x3C8E0A72F0539783),
        ), // 6.530612244897959183673469387755e-1
        Double::new(
            f64::from_bits(0x3FE4AFD6A052BF5A),
            f64::from_bits(0x3C90295FAD40A57F),
        ), // 6.464646464646464646464646464646e-1
        Double::new(
            f64::from_bits(0x3FE47AE147AE147A),
            f64::from_bits(0x3C9C28F5C28F5C29),
        ), // 6.400000000000000000000000000000e-1
        Double::new(
            f64::from_bits(0x3FE446F86562D9FA),
            f64::from_bits(0x3C9DC83CD4E93029),
        ), // 6.336633663366336633663366336634e-1
        Double::new(
            f64::from_bits(0x3FE4141414141414),
            f64::from_bits(0x3C64141414141414),
        ), // 6.274509803921568627450980392157e-1
        Double::new(
            f64::from_bits(0x3FE3E22CBCE4A902),
            f64::from_bits(0x3C8F1165E7254814),
        ), // 6.213592233009708737864077669903e-1
        Double::new(
            f64::from_bits(0x3FE3B13B13B13B13),
            f64::from_bits(0x3C96276276276276),
        ), // 6.153846153846153846153846153846e-1
        Double::new(
            f64::from_bits(0x3FE3813813813813),
            f64::from_bits(0x3C90270270270270),
        ), // 6.095238095238095238095238095238e-1
        Double::new(
            f64::from_bits(0x3FE3521CFB2B78C1),
            f64::from_bits(0x3C7A90E7D95BC60A),
        ), // 6.037735849056603773584905660377e-1
        Double::new(
            f64::from_bits(0x3FE323E34A2B10BF),
            f64::from_bits(0x3C89B8396BA9DE81),
        ), // 5.981308411214953271028037383178e-1
        Double::new(
            f64::from_bits(0x3FE2F684BDA12F68),
            f64::from_bits(0x3C82F684BDA12F68),
        ), // 5.925925925925925925925925925926e-1
        Double::new(
            f64::from_bits(0x3FE2C9FB4D812C9F),
            f64::from_bits(0x3C969B02593F69B0),
        ), // 5.871559633027522935779816513761e-1
        Double::new(
            f64::from_bits(0x3FE29E4129E4129E),
            f64::from_bits(0x3C804A7904A7904A),
        ), // 5.818181818181818181818181818182e-1
        Double::new(
            f64::from_bits(0x3FE27350B8812735),
            f64::from_bits(0x3C571024E6A17102),
        ), // 5.765765765765765765765765765766e-1
        Double::new(
            f64::from_bits(0x3FE2492492492492),
            f64::from_bits(0x3C82492492492492),
        ), // 5.714285714285714285714285714286e-1
        Double::new(
            f64::from_bits(0x3FE21FB78121FB78),
            f64::from_bits(0x3C621FB78121FB78),
        ), // 5.663716814159292035398230088496e-1
        Double::new(
            f64::from_bits(0x3FE1F7047DC11F70),
            f64::from_bits(0x3C81F7047DC11F70),
        ), // 5.614035087719298245614035087719e-1
        Double::new(
            f64::from_bits(0x3FE1CF06ADA2811C),
            f64::from_bits(0x3C9E0D5B450239E1),
        ), // 5.565217391304347826086956521739e-1
        Double::new(
            f64::from_bits(0x3FE1A7B9611A7B96),
            f64::from_bits(0x3C61A7B9611A7B96),
        ), // 5.517241379310344827586206896552e-1
        Double::new(
            f64::from_bits(0x3FE1811811811811),
            f64::from_bits(0x3C90230230230230),
        ), // 5.470085470085470085470085470085e-1
        Double::new(
            f64::from_bits(0x3FE15B1E5F75270D),
            f64::from_bits(0x3C415B1E5F75270D),
        ), // 5.423728813559322033898305084746e-1
        Double::new(
            f64::from_bits(0x3FE135C81135C811),
            f64::from_bits(0x3C7AE4089AE4089B),
        ), // 5.378151260504201680672268907563e-1
        Double::new(
            f64::from_bits(0x3FE1111111111111),
            f64::from_bits(0x3C61111111111111),
        ), // 5.333333333333333333333333333333e-1
        Double::new(
            f64::from_bits(0x3FE0ECF56BE69C8F),
            f64::from_bits(0x3C9BC4C2A50658DC),
        ), // 5.289256198347107438016528925620e-1
        Double::new(
            f64::from_bits(0x3FE0C9714FBCDA3A),
            f64::from_bits(0x3C982192E29F79B4),
        ), // 5.245901639344262295081967213115e-1
        Double::new(
            f64::from_bits(0x3FE0A6810A6810A6),
            f64::from_bits(0x3C90214D0214D021),
        ), // 5.203252032520325203252032520325e-1
        Double::new(
            f64::from_bits(0x3FE0842108421084),
            f64::from_bits(0x3C70842108421084),
        ), // 5.161290322580645161290322580645e-1
        Double::new(
            f64::from_bits(0x3FE0624DD2F1A9FB),
            f64::from_bits(0x3C9CED916872B021),
        ), // 5.120000000000000000000000000000e-1
        Double::new(
            f64::from_bits(0x3FE0410410410410),
            f64::from_bits(0x3C80410410410410),
        ), // 5.079365079365079365079365079365e-1
        Double::new(
            f64::from_bits(0x3FE0204081020408),
            f64::from_bits(0x3C60204081020408),
        ), // 5.039370078740157480314960629921e-1
    ];

    // x = 2^k * m, with 1 <= m < 2
    let k = x.hi().exponent();
    let m = x.hi().mant();

    // ln(m) = ln(hi) + ln(lo + 1)
    let hi = m >> (53 - 7) & 0x3F;
    let hi_f = f64::from_bits(0x3FF0_0000_0000_0000 | (hi << (53 - 7)));
    let hi = hi as usize;
    let scale = scalbn(1.0, (-k).into());
    let lo = (x.pmul1(scale) - hi_f) * LN_LO_SCALE_TBL[hi];
    let ln_hi = LN_TBL[hi];

    // The upper half of `LN_TBL` holds `ln((1 + i / 64) / 2)` instead of
    // `ln(1 + i / 64)` to avoid cancellation when `x` is close to 1.
    let k = if hi >= 32 { k + 1 } else { k };

    (f64::from(k), lo, ln_hi)
}
