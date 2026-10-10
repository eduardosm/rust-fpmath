//! Logarithmic functions for `f64`, and the parts that `f32` shares.
//!
//! The argument is reduced to `x = 2^k * (1 + z) / s`, where `s ~= 1 / c`,
//! with `c = 1 + i / N` and `N = 2^LN_TBL_BITS`, is taken from `LN_TBL` and
//! `|z| <= 1 / (2 * N)`, so
//!
//! `ln(x) = (k * ln(2) - ln(s)) + ln(1 + z)`
//!
//! `s` has 24 significant bits, so `z` is calculated exactly as a sum of two
//! `f64`, `zh + zl`, splitting the mantissa of `x` in two halves.
//!
//! `-ln(s)` is also taken from `LN_TBL`, as `hi + lo`, where `hi` is a
//! multiple of 2^-42, so `t_hi = k * ln(2)_hi - ln(s)_hi` is exact
//! (`ln(2)_hi + ln(2)_lo` is the entry of `LN_TBL[N]`, where `s = 1/2`).
//! `t_hi` is zero when `x` is close to one (`(k, i)` is `(0, 0)` or
//! `(-1, N)`), and its magnitude is at least about twice the one of
//! `ln(1 + z)` otherwise, so `t_hi + ln(1 + z)` has no cancellation.
//!
//! `ln(1 + z)` is approximated with a polynomial, where `z - z^2 / 2` is
//! calculated with extended precision, so the relative error is small even
//! when `t_hi` is zero (see `ln_1p_poly`).
//!
//! `log2(x)` and `log10(x)` are calculated as `ln(x) * log2(e)` and
//! `ln(x) * log10(e)`, with extended precision.

use super::f64x2::F64x2;
use crate::traits::Float as _;

/// Number of bits of the index of `LN_TBL`
pub(crate) const LN_TBL_BITS: u32 = {
    let n = LN_TBL.len() - 1;
    assert!(n.is_power_of_two());
    n.ilog2()
};

/// `N`
const TBL_N: usize = 1 << LN_TBL_BITS;

// `log2(e)` and `log10(e)`, the first parts have 26 bits, so their products
// with numbers of up to 27 bits are exact
// GENERATE: split_const LOG2_E_X LOG2_E 0 26 53
const LOG2_E_X0: f64 = f64::from_bits(0x3FF7154760000000); // 1.4426950216293335e0
const LOG2_E_X1: f64 = f64::from_bits(0x3E54AE0BF85DDF44); // 1.9259629911266175e-8

// GENERATE: split_const LOG10_E_X LOG10_E 0 26 53
const LOG10_E_X0: f64 = f64::from_bits(0x3FDBCB7B10000000); // 4.342944771051407e-1
const LOG10_E_X1: f64 = f64::from_bits(0x3E349B9438CA9AAE); // 4.798111141615973e-9

// GENERATE: ln_table f32 7 42
// LN_TBL[i] = (bits(s), bits(hi), bits(lo)), where s = 1 / (1 + i / 128) rounded to f32
// and hi + lo = -ln(s), with hi being a multiple of 2^-42
static LN_TBL: [(u64, u64, u64); 129] = [
    (0x3FF0000000000000, 0x8000000000000000, 0x0000000000000000), // -0
    (0x3FEFC07F00000000, 0x3F7FE02B6B100000, 0x3D19E43F0DDA563A), // 7.782144167345254348270883701295e-3
    (0x3FEF81F820000000, 0x3F8FC0A890FC0000, 0x3CDF207CF6D3A147), // 1.550418560464267996905639914238e-2
    (0x3FEF4465A0000000, 0x3F97B91ACFD60000, 0xBD33B8F3B602B076), // 2.316705602190537238721498990905e-2
    (0x3FEF07C200000000, 0x3F9F82990E780000, 0x3D29C0267C68B48F), // 3.077162886443174476491673439972e-2
    (0x3FEECC07C0000000, 0x3FA39E86E1FE8000, 0x3D3EC69C80A727D5), // 3.831883915642740072918360008840e-2
    (0x3FEE9131A0000000, 0x3FA77459BE330000, 0xBD316E54E58198F4), // 4.580955931435883960418948687401e-2
    (0x3FEE573AC0000000, 0x3FAB42DE09198000, 0xBD1C555AE5CD81F7), // 5.324453221394135711875754445375e-2
    (0x3FEE1E1E20000000, 0x3FAF0A30A0118000, 0xBD2D589E8336993C), // 6.062461809114455105758595621472e-2
    (0x3FEDE5D6E0000000, 0x3FB1653710A38000, 0xBD147356768ED653), // 6.795066982474966495957561077828e-2
    (0x3FEDAE6080000000, 0x3FB341D7461BC000, 0x3D31DD129980DB66), // 7.522340261113620686138047104822e-2
    (0x3FED77B660000000, 0x3FB51B06DD060000, 0x3D38522B27899EE8), // 8.244364639367177350587408978777e-2
    (0x3FED41D420000000, 0x3FB6F0D272E58000, 0xBD34B3441B665813), // 8.961215310175170053959160125861e-2
    (0x3FED0CB580000000, 0x3FB8C3465E318000, 0x3D3B4515ACC0F5BB), // 9.672965812351915055693572820775e-2
    (0x3FECD85680000000, 0x3FBA926D8A4AC000, 0x3D356FE50BD4C547), // 1.037968123080952306079818680609e-1
    (0x3FECA4B300000000, 0x3FBC5E54BF5BC000, 0x3D1D1E575861FE06), // 1.108143775161610720305939572774e-1
    (0x3FEC71C720000000, 0x3FBE27074E2B0000, 0xBD2A302C2AF05591), // 1.177830282058028853705414622360e-1
    (0x3FEC3F8F00000000, 0x3FBFEC9141DC0000, 0xBD3544D5D1AE60B1), // 1.247034822262475412642144990210e-1
    (0x3FEC0E0700000000, 0x3FC0D77E8CD08000, 0x3D3CB4CD2EE31F2C), // 1.315763652392998972681200071818e-1
    (0x3FEBDD2B80000000, 0x3FC1B72B012F6000, 0x3D2E9EE418189241), // 1.384023433482159871273966776740e-1
    (0x3FEBACF920000000, 0x3FC29552C4200000, 0xBD35A447F44CD6A7), // 1.451819856301112504477563649778e-1
    (0x3FEB7D6C40000000, 0x3FC371FC161E8000, 0x3D3EE93F9B2D8052), // 1.519160373692291128364325379689e-1
    (0x3FEB4E81C0000000, 0x3FC44D2B38CB8000, 0xBD16B841614C5AE7), // 1.586050059622519372591924758631e-1
    (0x3FEB203640000000, 0x3FC526E5E5A1C000, 0xBD3790B237FC5223), // 1.652495738266297379247708341657e-1
    (0x3FEAF286C0000000, 0x3FC5FF3060A7A000, 0xBD38566F183C169C), // 1.718502494760786531718462988206e-1
    (0x3FEAC57020000000, 0x3FC6D60FCE19E000, 0xBD3BC2035713EA29), // 1.784076458312861821885540506352e-1
    (0x3FEA98EF60000000, 0x3FC7AB890410E000, 0xBD2BDB8072534A2D), // 1.849223394253345677130629775477e-1
    (0x3FEA6D01A0000000, 0x3FC87FA08620C000, 0x3D3229A240137954), // 1.913948683664520538343121086151e-1
    (0x3FEA41A420000000, 0x3FC9525A80F46000, 0xBD3290F37D9FFA39), // 1.978257302914039207473224364487e-1
    (0x3FEA16D400000000, 0x3FCA23BBFFE2C000, 0xBD3531CD91DDF460), // 2.042155265275298086784654958044e-1
    (0x3FE9EC8EA0000000, 0x3FCAF3C91880C000, 0xBC8C331A31AE8320), // 2.105647439616404392049810876970e-1
    (0x3FE9C2D140000000, 0x3FCBC286BE2D8000, 0x3D39D71BF3AD8F32), // 2.168739727595502141009178959878e-1
    (0x3FE99999A0000000, 0x3FCC8FF7A79AA000, 0xBD27694F68A22EDF), // 2.231435364130486729409401999180e-1
    (0x3FE970E500000000, 0x3FCD5C21434FC000, 0xBD21A191BBCF9D71), // 2.293740824383945111542558561107e-1
    (0x3FE948B100000000, 0x3FCE27075E2B0000, 0xBD3A322C2AF02AE7), // 2.355660638621863399093355717065e-1
    (0x3FE920FB40000000, 0x3FCEF0ADFDDC6000, 0xBD2AFFA79C7C82F9), // 2.417199601702098045818177351729e-1
    (0x3FE8F9C180000000, 0x3FCFB918BD5E4000, 0xBD0BC72AAAF291DC), // 2.478362011574849352891510143993e-1
    (0x3FE8D30180000000, 0x3FD04025B6B4D000, 0x3CF278B89FC0E2D5), // 2.539152416459314823986879696269e-1
    (0x3FE8ACB900000000, 0x3FD0A3250A739000, 0x3D0DFBEE7F9AADB9), // 2.599575616898297454806277431951e-1
    (0x3FE886E600000000, 0x3FD1058BD1AE5000, 0xBD34799D81922822), // 2.659635112442356506093584593314e-1
    (0x3FE8618620000000, 0x3FD1675C97ABA000, 0x3D38448E731CBB19), // 2.719336968571904399944446255976e-1
    (0x3FE83C9780000000, 0x3FD1C898B369A000, 0xBD380DF0E5C70FAA), // 2.778684379649403465711003981751e-1
    (0x3FE8181820000000, 0x3FD22941E6CF8000, 0xBD3A5BAEF5EE0D23), // 2.837681535728707226751131270095e-1
    (0x3FE7F40600000000, 0x3FD2895A0BDE8000, 0x3D3A8F7AD24BE946), // 2.896332851324621077106404082911e-1
    (0x3FE7D05F40000000, 0x3FD2E8E2BEE12000, 0xBD267A1E99B7212D), // 2.954642166191261817874898896951e-1
    (0x3FE7AD2200000000, 0x3FD347DDB2988000, 0xBD25354DD4BC8092), // 3.012613529299038215845441762068e-1
    (0x3FE78A4C80000000, 0x3FD3A64C59694000, 0x3D37A79CBCD73B26), // 3.070250390202021674759325241755e-1
    (0x3FE767DCE0000000, 0x3FD404309206A000, 0x3D3F9316304A7690), // 3.127557207141065538185448233767e-1
    (0x3FE745D180000000, 0x3FD4618BA21C6000, 0xBD13582F48772F77), // 3.184537013162126722041357403935e-1
    (0x3FE7242880000000, 0x3FD4BE5F93778000, 0xBD3D7C72CD9AD8CF), // 3.241194667915668285944372030077e-1
    (0x3FE702E060000000, 0x3FD51AAD7C2E0000, 0xBD3F4810DB0AEBAC), // 3.297532761279197135195439299531e-1
    (0x3FE6E1F760000000, 0x3FD5767736C56000, 0xBD362FAB951AAB22), // 3.353555712577993609646034975968e-1
    (0x3FE6C16C20000000, 0x3FD5D1BDA5581000, 0xBD38C19DC9CD7AE3), // 3.409265627562065634709105010176e-1
    (0x3FE6A13CE0000000, 0x3FD62C82C939C000, 0x3D3E8A8FBD65467B), // 3.464667286963224812829611700087e-1
    (0x3FE6816820000000, 0x3FD686C8039B1000, 0x3D32D1D90AF1D814), // 3.519763980114689861908757312240e-1
    (0x3FE661EC60000000, 0x3FD6E08EC7ABA000, 0x3D1EA5893952FAE7), // 3.574559163958201027934088559424e-1
    (0x3FE642C860000000, 0x3FD739D7E2BBD000, 0x3CC379C4975AA053), // 3.629054750629171343005994770421e-1
    (0x3FE623FA80000000, 0x3FD792A545DD4000, 0x3D3E9F105763673F), // 3.683255369443210062140514552876e-1
    (0x3FE6058160000000, 0x3FD7EAF83C82B000, 0xBCEE4CA62D0C2303), // 3.737164107249066558701762176043e-1
    (0x3FE5E75BC0000000, 0x3FD842D1C51E9000, 0xBD33951313B16C3C), // 3.790783333771955827190652502903e-1
    (0x3FE5C98820000000, 0x3FD89A33A8C14000, 0x3D231DED38DD9F2D), // 3.844117305753000779961547342016e-1
    (0x3FE5AC0560000000, 0x3FD8F11EA7B66000, 0x3D267DF4BB6504A7), // 3.897167814080093464489424802880e-1
    (0x3FE58ED240000000, 0x3FD94793EE211000, 0x3D3C2D093354A29F), // 3.949937654000314634630749002152e-1
    (0x3FE571ED40000000, 0x3FD99D957617E000, 0x3D0177B525DA119B), // 4.002431538824644386344462200275e-1
    (0x3FE5555560000000, 0x3FD9F323CCBFA000, 0xBD3EB03525D4D0EE), // 4.054650783058424383719016422673e-1
    (0x3FE5390940000000, 0x3FDA4840ABE5C000, 0xBD33C0880D5EBB85), // 4.106599501309782167055849217576e-1
    (0x3FE51D07E0000000, 0x3FDA9CECBB9A1000, 0xBD3EB7C2DEB47883), // 4.158279258773564002025959198250e-1
    (0x3FE5015020000000, 0x3FDAF12910C78000, 0xBD3E30A931FB149B), // 4.209692634448238732086819937118e-1
    (0x3FE4E5E0A0000000, 0x3FDB44F791CC9000, 0xBD033568222EE824), // 4.260844157999969145666156836678e-1
    (0x3FE4CAB880000000, 0x3FDB9858AC931000, 0x3D0FE431F645ABC9), // 4.311734853074681923012432830279e-1
    (0x3FE4AFD6A0000000, 0x3FDBEB4D9EA72000, 0xBD321019E78B213C), // 4.362367677062406453982007079496e-1
    (0x3FE49539E0000000, 0x3FDC3DD7B34DB000, 0xBD258C1E61F4A6B1), // 4.412745715150848949217947053497e-1
    (0x3FE47AE140000000, 0x3FDC8FF7DF9AA000, 0xBD37674F689B0434), // 4.462871249801615521042588185927e-1
    (0x3FE460CBC0000000, 0x3FDCE1AF2485F000, 0x3D2F82EBE9688193), // 4.512746674225232215822026023384e-1
    (0x3FE446F860000000, 0x3FDD32FE8F00F000, 0xBD30A17084DB36E6), // 4.562374493140714881777127660099e-1
    (0x3FE42D6620000000, 0x3FDD83E7380A3000, 0xBD07E065D47B2558), // 4.611757323516379451816315828079e-1
    (0x3FE4141420000000, 0x3FDDD469DEC1C000, 0x3D32B01B9888B5CA), // 4.660896945343420154055957062997e-1
    (0x3FE3FB0140000000, 0x3FDE2488197C7000, 0xBD2ECF0A1385D380), // 4.709797142874684383651002091563e-1
    (0x3FE3E22CC0000000, 0x3FDE744257D68000, 0x3D3E22ADF68D699E), // 4.758458955567382114785110602100e-1
    (0x3FE3C995A0000000, 0x3FDEC399E0C69000, 0xBD29F221188B644B), // 4.806885428499293307824604589280e-1
    (0x3FE3B13B20000000, 0x3FDF128F37AF0000, 0x3D3BE4CD71F9EEF7), // 4.855077785287985170720236099847e-1
    (0x3FE3991C20000000, 0x3FDF612421F03000, 0xBD3D1DB8EA3AFD52), // 4.903040248324362121123399697759e-1
    (0x3FE3813820000000, 0x3FDFAF586678F000, 0x3D295FDD7D72487F), // 4.950772286136266843808676654914e-1
    (0x3FE3698E00000000, 0x3FDFFD2DE057F000, 0x3D3293565F2C03DD), // 4.998278323035470390915639480409e-1
    (0x3FE3521D00000000, 0x3FE025529DA5D000, 0x3D1FF8D38D265A88), // 5.045559958512342042329536413464e-1
    (0x3FE33AE460000000, 0x3FE04BDF95E92800, 0xBD22CA9C33F263F3), // 5.092618873543081444559857001904e-1
    (0x3FE323E340000000, 0x3FE0723E6D1CE000, 0xBD1765B50D05C088), // 5.139457827672023550623707599945e-1
    (0x3FE30D1900000000, 0x3FE0986F51573800, 0xBD36F9B7012B52B1), // 5.186077679333359375537849800050e-1
    (0x3FE2F684C0000000, 0x3FE0BE72E0252800, 0x3D3417B4C4BDAEF4), // 5.232481363139672673485545777003e-1
    (0x3FE2E025C0000000, 0x3FE0E4498651D000, 0xBD3BA040A8D10B36), // 5.278670905521649601630516026668e-1
    (0x3FE2C9FB40000000, 0x3FE109F3B52D5000, 0xBD3B05F0E07B784F), // 5.324648417103111938546804430735e-1
    (0x3FE2B404A0000000, 0x3FE12F71ABD3F000, 0xBCFDF85F6CA82541), // 5.370415073407390837520877431278e-1
    (0x3FE29E4120000000, 0x3FE154C3E3F4D800, 0xBD308E93865617F8), // 5.415973140977124098379069756207e-1
    (0x3FE288B020000000, 0x3FE179EAA4989800, 0x3D2A932060A58498), // 5.461323942916368684890781422820e-1
    (0x3FE27350C0000000, 0x3FE19EE6A767C800, 0x3D271705CBD2062D), // 5.506470937382756324257694804422e-1
    (0x3FE25E2280000000, 0x3FE1C3B804714000, 0xBD3E7EFB586BDB02), // 5.551414572490838280930121033911e-1
    (0x3FE24924A0000000, 0x3FE1E85F46704000, 0x3D1B27BD8AA4BE7D), // 5.596157432319401039286121346667e-1
    (0x3FE2345680000000, 0x3FE20CDCC492A800, 0x3D3B81BA81E2C303), // 5.640701140704163192372055913274e-1
    (0x3FE21FB780000000, 0x3FE23130D9BEC000, 0xBD17ADA4392F0651), // 5.685047390779590174795467485071e-1
    (0x3FE20B4700000000, 0x3FE2555BE498F800, 0xBCF699FDE0D6ECD3), // 5.729197945399796317799967239307e-1
    (0x3FE1F70480000000, 0x3FE2795E0E89B000, 0x3D21B2B783F38641), // 5.773153575842430351498594142850e-1
    (0x3FE1E2EF40000000, 0x3FE29D37F642B000, 0x3D117D2B9AD30F0F), // 5.816917238021388393912458045281e-1
    (0x3FE1CF06A0000000, 0x3FE2C0EA05C49000, 0xBD26AE5855F04C33), // 5.860490906383854063303648356937e-1
    (0x3FE1BB4A40000000, 0x3FE2E47437640000, 0x3D33420AA10C34A6), // 5.903874475334989496910760930121e-1
    (0x3FE1A7B960000000, 0x3FE307D7354F1000, 0x3D17C5F6B2145402), // 5.947071114719830949151515301694e-1
    (0x3FE1945380000000, 0x3FE32B133A122000, 0xBD34764FD54A4B7C), // 5.990081915087285503472805338426e-1
    (0x3FE1811820000000, 0x3FE34E28831CE000, 0x3D2DE0FEB8CC9B88), // 6.032908020779890259297447759563e-1
    (0x3FE16E0680000000, 0x3FE37117C6474800, 0xBD016D88BF07941E), // 6.075552818895098337624498245114e-1
    (0x3FE15B1E60000000, 0x3FE393E0D2562800, 0x3D30CD6E2213010C), // 6.118015392433477560336562090015e-1
    (0x3FE1485F00000000, 0x3FE3B68464000000, 0xBD3E960388DC2E7E), // 6.160299256442890723255563617656e-1
    (0x3FE135C820000000, 0x3FE3D9024EF15800, 0xBD1F315F7C1100FA), // 6.202403585291172368847601817694e-1
    (0x3FE12358E0000000, 0x3FE3FB5B92917000, 0xBD17560E2C3AE020), // 6.244333136232646309393607004680e-1
    (0x3FE1111120000000, 0x3FE41D8FCC467000, 0x3D35D52325ACECF8), // 6.286086072683113193006692094344e-1
    (0x3FE0FEF020000000, 0x3FE43F9FC6B9D000, 0xBD28C0D4D8BBB64B), // 6.327666169513137481889558173122e-1
    (0x3FE0ECF560000000, 0x3FE4618BD89C6000, 0xBD23599F4811314C), // 6.369075041465859675208118802952e-1
    (0x3FE0DB20A0000000, 0x3FE48353E22A8800, 0x3D1C7282BD4418B9), // 6.410312096889154241340802106899e-1
    (0x3FE0C97140000000, 0x3FE4A4F87BB04000, 0xBD236EDD837EE591), // 6.451380172529407398451259739004e-1
    (0x3FE0B7E6E0000000, 0x3FE4C679C70CF000, 0xBD2BDDC570F0B777), // 6.492279899316104762370770022819e-1
    (0x3FE0A68100000000, 0x3FE4E7D825B75800, 0x3D3DB6E584D78782), // 6.533013092656493172671641298636e-1
    (0x3FE0953F40000000, 0x3FE50913BE816800, 0x3D0B8932CE6380D6), // 6.573580475626508316773182982105e-1
    (0x3FE0842100000000, 0x3FE52A2D365BC800, 0xBD32888C41AFDCA8), // 6.613985120476878400447670120323e-1
    (0x3FE0732600000000, 0x3FE54B247B999800, 0xBD3B10B755D6D08C), // 6.654226697979941274586408586594e-1
    (0x3FE0624DE0000000, 0x3FE56BF9BC33F000, 0x3D3D1B50E2E58B72), // 6.694306064451790899133855097664e-1
    (0x3FE0519800000000, 0x3FE58CADA5CD7800, 0x3D28D3092F1083DB), // 6.734226454098447766918489156040e-1
    (0x3FE0410420000000, 0x3FE5AD402D35A000, 0xBD18801ACBE194A5), // 6.773987658498081817192108816867e-1
    (0x3FE03091C0000000, 0x3FE5CDB1C6EC1800, 0xBD12874123EECB74), // 6.813591847610331623584007495584e-1
    (0x3FE0204080000000, 0x3FE5EE02AB241800, 0xBD28A7F29F69F831), // 6.853040068242097219448560626074e-1
    (0x3FE0101020000000, 0x3FE60E32D4878800, 0x3D1D11578FDDDC20), // 6.892332220998272409516436528915e-1
    (0x3FE0000000000000, 0x3FE62E42FEFA3800, 0x3D2EF35793C76730), // 6.931471805599453094172321214582e-1
];

impl crate::generic::Log for f64 {
    #[inline]
    fn ln_finite(x: Self, edelta: i16) -> Self {
        let (hi, lo) = ln_parts(x, edelta);
        hi + lo
    }

    #[inline]
    fn ln_1p_finite(x: Self) -> Self {
        let x_exp = x.exponent();
        if x_exp < -8 {
            if x_exp < -53 {
                // ln(1 + x) = x - x^2 / 2 + ..., where |x| < 2^-53, so
                // x^2 / 2 - ... < 2^-54 * |x|, which is less than half the
                // spacing of `f64` around x, and the result rounds to x
                return x;
            }
            // |x| < 2^-8
            let (hi, lo) = ln_1p_poly(x, 0.0);
            hi + lo
        } else if x_exp >= 64 {
            // ln(1 + x) = ln(x) + ln(1 + 1 / x), where ln(1 + 1 / x) < 2^-64
            // is less than 2^-69 relative to ln(x) > 44, so it is negligible
            // (and calculating it could involve subnormal numbers)
            let (hi, lo) = ln_parts(x, 0);
            hi + lo
        } else {
            // 1 + x = p_hi + p_lo
            let p = F64x2::add11(1.0, x);
            let (hi, lo) = ln_sum_parts(p.hi(), p.lo());
            hi + lo
        }
    }

    #[inline]
    fn log2_finite(x: Self, edelta: i16) -> Self {
        let (hi, lo) = ln_parts(x, edelta);
        mul_split(LOG2_E_X0, LOG2_E_X1, hi, lo)
    }

    #[inline]
    fn log10_finite(x: Self, edelta: i16) -> Self {
        let (hi, lo) = ln_parts(x, edelta);
        mul_split(LOG10_E_X0, LOG10_E_X1, hi, lo)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(x * 2^edelta)`, with the
/// relative error of `ln_1p_poly` and `|lo| <= ~2^-15 * |hi|`.
///
/// `x` must be positive and normal.
#[inline]
pub(super) fn ln_parts(x: f64, edelta: i16) -> (f64, f64) {
    let red = Reduced::new(x, edelta);
    let (hi, lo) = ln_1p_poly(red.zh, red.zl);
    red.add_t(hi, lo)
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(t_hi + t_lo)`, with the
/// relative error of `ln_1p_poly` (about 2^-66) when
/// `|ln(t_hi + t_lo)| > 2^-37`.
///
/// `ln(1 + t_lo / t_hi)` is approximated by `t_lo / t_hi`, with an absolute
/// error of up to 2^-101, which is not negligible for smaller results (for
/// example, the result is not zero when `t_hi + t_lo = 1` with `t_lo != 0`).
///
/// `t_hi` must be positive and normal, and `|t_lo| <= 2^-50 * t_hi`.
#[inline]
pub(super) fn ln_sum_parts(t_hi: f64, t_lo: f64) -> (f64, f64) {
    // ln(t_hi + t_lo) = ln(t_hi) + ln(1 + d), where d = t_lo / t_hi and
    // ln(1 + d) = d - d^2 / 2 + ..., with |d^2 / 2| <= 2^-101
    let (hi, lo) = ln_parts(t_hi, 0);
    (hi, lo + t_lo / t_hi)
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(x * 2^edelta)`, with a
/// relative error of about 2^-75.
///
/// It is slower than `ln_parts`, and meant for the functions that amplify the
/// error of the logarithm: `pow`, by up to `|y * ln(x)|` (about 2^9.5), and
/// the gamma functions (see `super::gamma`).
///
/// `x` must be positive and normal.
#[inline]
pub(super) fn ln_accurate_parts(x: f64, edelta: i16) -> (f64, f64) {
    let red = Reduced::new(x, edelta);
    let (hi, lo) = ln_1p_poly_accurate(red.zh, red.zl);
    red.add_t(hi, lo)
}

/// Like `ln_sum_parts` (including its absolute error), but with the relative
/// error of `ln_accurate_parts`.
#[inline]
pub(super) fn ln_accurate_sum_parts(t_hi: f64, t_lo: f64) -> (f64, f64) {
    let (hi, lo) = ln_accurate_parts(t_hi, 0);
    (hi, lo + t_lo / t_hi)
}

/// Splits `x * 2^edelta = 2^k * m`, with `1 <= m < 2`, returning `(k, m, i)`,
/// where `i = round((m - 1) * N)`, so `0 <= i <= N`.
///
/// `x` must be positive and normal.
#[inline]
pub(crate) fn split_ln_arg(x: f64, edelta: i16) -> (i16, f64, usize) {
    let k = x.exponent() + edelta;
    let bits = x.to_bits();
    let m = f64::from_bits((bits & f64::MANT_MASK) | 1.0f64.to_bits());
    // round to nearest (ties up) with one more bit
    let i2 = (bits >> (f64::MANT_BITS - LN_TBL_BITS as u8 - 1)) as usize & ((2 * TBL_N) - 1);
    (k, m, (i2 + 1) >> 1)
}

/// Returns `(s, t_hi, t_lo)`, where `s ~= 1 / (1 + i / N)` (with 24 bits) and
/// `t_hi + t_lo ~= k * ln(2) - ln(s)`, with a relative error of at most
/// 2^-85.
///
/// `t_hi` is exact (`|k| < 2^11`), and it is zero when `(k, i)` is `(0, 0)`
/// or `(-1, N)`.
#[inline]
pub(crate) fn ln_tbl(k: i16, i: usize) -> (f64, f64, f64) {
    let k = f64::from(k);
    let (s, tbl_hi, tbl_lo) = LN_TBL[i];
    // `s = 1/2` in `LN_TBL[N]`, so it holds `-ln(1/2) = ln(2)`
    let (_, ln2_hi, ln2_lo) = LN_TBL[TBL_N];
    // exact: `ln2_hi` and `tbl_hi` are multiples of 2^-42
    let t_hi = k * f64::from_bits(ln2_hi) + f64::from_bits(tbl_hi);
    let t_lo = k * f64::from_bits(ln2_lo) + f64::from_bits(tbl_lo);
    (f64::from_bits(s), t_hi, t_lo)
}

/// Reduced argument: `x = 2^k * (1 + zh + zl) / s`
struct Reduced {
    t_hi: f64,
    t_lo: f64,
    zh: f64,
    zl: f64,
}

impl Reduced {
    /// Reduces `x * 2^edelta`, where `x` must be positive and normal.
    ///
    /// `zh + zl` is exact, `|zh + zl| <= 2^-8` and `|zl| <= ulp(zh) / 2`.
    #[inline]
    fn new(x: f64, edelta: i16) -> Self {
        let (k, m, i) = split_ln_arg(x, edelta);
        let (s, t_hi, t_lo) = ln_tbl(k, i);

        // m = mh + ml, where `mh` has 26 bits and `ml` at most 27
        let mh = m.split_hi();
        let ml = m - mh;
        // `mh * s` and `ml * s` are exact (`s` has 24 bits), and so is
        // `mh * s - 1` (by Sterbenz lemma).
        let z1 = mh * s - 1.0;
        let z2 = ml * s;
        // `z1` is a multiple of 2^-49 and `z2` is a multiple of 2^-76, with
        // `|z2| < 2^-25`. If `|z1 + z2| < 2^-24`, the sum is exact (it fits in
        // 52 bits), otherwise `|z1| > 2^-25 > |z2|`, so the sum is exact as
        // a sum of two `f64`.
        let z = F64x2::fast_add11(z1, z2);

        Self {
            t_hi,
            t_lo,
            zh: z.hi(),
            zl: z.lo(),
        }
    }

    /// Returns `(hi, lo)` such that `hi + lo ~= k * ln(2) - ln(s) + p_hi + p_lo`,
    /// where `p_hi + p_lo ~= ln(1 + zh + zl)`.
    #[inline]
    fn add_t(&self, p_hi: f64, p_lo: f64) -> (f64, f64) {
        // `t_hi` is zero or larger than `p_hi` in magnitude
        let s = F64x2::fast_add11(self.t_hi, p_hi);
        (s.hi(), (s.lo() + self.t_lo) + p_lo)
    }
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(1 + zh + zl)`, for
/// `|zh + zl| <= 2^-8` and `|zl| <= 2^-51 * |zh|`, with a relative error of
/// about 2^-66:
/// * Polynomial approximation: 2^-71.4 (see `ln_1p_q`).
/// * Evaluation of the terms after `zh - zh^2 / 2` (magnitude below
///   2^-16 * |zh|): ~2^-67.
/// * Neglected `zl` terms: 2^-67.
///
/// `|lo| <= 2^-16 * |hi|`.
#[inline]
fn ln_1p_poly(zh: f64, zl: f64) -> (f64, f64) {
    // ln(1 + zh + zl) ~= ln(1 + zh) + zl * (1 - zh)
    // ln(1 + zh) ~= zh - zh^2 / 2 + zh^3 * Q(zh)
    //
    // zh = a + b, where `a` has 26 bits, so `a^2 / 2` is exact
    // zh^2 / 2 = a^2 / 2 + b * (zh + a) / 2
    let a = zh.split_hi();
    let b = zh - a;
    let s = F64x2::fast_add11(zh, -(0.5 * a) * a);

    let zh2 = zh * zh;
    let lo = (s.lo() + (zl - zl * zh)) - (0.5 * b) * (zh + a);
    (s.hi(), lo + (zh2 * zh) * ln_1p_q(zh))
}

/// Returns `Q(z)` such that `ln(1 + z) ~= z - z^2 / 2 + z^3 * Q(z)`, for
/// `|z| <= 2^-8`, with a relative error of 2^-71.4 (in `ln(1 + z)`).
#[inline]
pub(crate) fn ln_1p_q(z: f64) -> f64 {
    // GENERATE: ln_1p_poly f64 6 -0.003907 0.003907 2
    const K3: f64 = f64::from_bits(0x3FD5555555555557); // 3.333333333333334e-1
    const K4: f64 = f64::from_bits(0xBFD0000000000003); // -2.5000000000000017e-1
    const K5: f64 = f64::from_bits(0x3FC999999989C304); // 1.9999999997119045e-1
    const K6: f64 = f64::from_bits(0xBFC55555553F1CA4); // -1.6666666662624607e-1
    const K7: f64 = f64::from_bits(0x3FC2493FE50367B0); // 1.428604000309952e-1
    const K8: f64 = f64::from_bits(0xBFC0001DCEA16A3B); // -1.2500355328923426e-1

    let z2 = z * z;
    (K3 + z * K4) + z2 * ((K5 + z * K6) + z2 * (K7 + z * K8))
}

/// Returns `(hi, lo)` such that `hi + lo ~= ln(1 + zh + zl)`, for
/// `|zh + zl| <= 2^-8` and `|zl| <= 2^-52 * |zh|`, with a relative error of
/// about 2^-77:
/// * Polynomial approximation: 2^-79.4.
/// * Evaluation of the terms after `zh - zh^2 / 2 + zh^3 / 3` (magnitude below
///   2^-25.9 * |zh|): ~2^-77.
/// * Neglected `zl` terms: 2^-77.
///
/// `|lo| <= 2^-25 * |hi|`.
#[inline]
fn ln_1p_poly_accurate(zh: f64, zl: f64) -> (f64, f64) {
    // ln(1 + z) - z + z^2 / 2 - z^3 / 3 ~= z^4 * (K4 + K5 * z + ... + K9 * z^5)
    // GENERATE: ln_1p_poly f64 6 -0.003907 0.003907 3
    const K4: f64 = f64::from_bits(0xBFD0000000000002); // -2.500000000000001e-1
    const K5: f64 = f64::from_bits(0x3FC999999999999F); // 2.0000000000000015e-1
    const K6: f64 = f64::from_bits(0xBFC5555555449C4C); // -1.6666666663624807e-1
    const K7: f64 = f64::from_bits(0x3FC2492492330E38); // 1.4285714281696626e-1
    const K8: f64 = f64::from_bits(0xBFC0001A3E69B0CA); // -1.2500312850477818e-1
    const K9: f64 = f64::from_bits(0x3FBC71FF703FF33D); // 1.1111446848367464e-1

    // 1/3 = K3_HI + K3_LO, where `K3_HI` has 11 bits and `K3_LO` has a
    // relative error of 2^-53 (2^-68.6 relative to 1/3)
    const K3_HI: f64 = f64::from_bits((1.0f64 / 3.0).to_bits() & (u64::MAX << 42));
    const K3_LO: f64 = (1.0 - 3.0 * K3_HI) / 3.0;

    // ln(1 + zh + zl) ~= ln(1 + zh) + zl * (1 - zh + zh^2)
    // ln(1 + zh) ~= zh - zh^2 / 2 + zh^3 / 3 + zh^4 * Q(zh)
    //
    // The leading parts of the quadratic and cubic terms are added exactly to
    // `zh`, so the rounding errors of the remaining terms are small:
    // * zh = a2 + b2, where `a2` has 26 bits, so `a2^2 / 2` is exact, and
    //   zh^2 / 2 = a2^2 / 2 + b2 * (zh + a2) / 2
    // * zh = a3 + b3, where `a3` has 14 bits, so `K3_HI * a3^3` is exact, and
    //   zh^3 / 3 = K3_HI * a3^3 + K3_HI * b3 * (zh^2 + zh * a3 + a3^2) + K3_LO * zh^3
    let a2 = zh.split_hi();
    let b2 = zh - a2;
    let a3 = f64::from_bits(zh.to_bits() & (u64::MAX << 39));
    let b3 = zh - a3;
    let s1 = F64x2::fast_add11(zh, -(0.5 * a2) * a2);
    let s2 = F64x2::fast_add11(s1.hi(), K3_HI * ((a3 * a3) * a3));

    let zh2 = zh * zh;
    let zh4 = zh2 * zh2;
    let q = ((K4 + zh * K5) + zh2 * (K6 + zh * K7)) + zh4 * (K8 + zh * K9);
    let c2 = (0.5 * b2) * (zh + a2);
    let c3 = K3_HI * (b3 * ((zh2 + zh * a3) + a3 * a3)) + K3_LO * (zh2 * zh);
    let cl = zl * ((1.0 - zh) + zh2);
    let lo = ((((s1.lo() + s2.lo()) - c2) + c3) + cl) + zh4 * q;
    (s2.hi(), lo)
}

/// Returns `(c_hi + c_lo) * (u_hi + u_lo)` rounded to `f64`, where `c_hi` has
/// at most 26 significant bits, `|c_lo| <= 2^-26 * |c_hi|` and
/// `|u_lo| <= ~2^-15 * |u_hi|`, with a relative error of about 2^-66 before
/// the final rounding.
#[inline]
fn mul_split(c_hi: f64, c_lo: f64, u_hi: f64, u_lo: f64) -> f64 {
    // u_hi = u1 + u2, where `u1` has 26 bits and `u2` at most 27, so
    // `c_hi * u1` and `c_hi * u2` are exact
    // (c_hi + c_lo) * (u_hi + u_lo) = c_hi * u1 + c_hi * u2 + c_hi * u_lo + c_lo * (u_hi + u_lo)
    let u1 = u_hi.split_hi();
    let u2 = u_hi - u1;
    let hi = c_hi * u1;
    let lo = c_hi * u2 + (c_lo * (u_hi + u_lo) + c_hi * u_lo);
    hi + lo
}
