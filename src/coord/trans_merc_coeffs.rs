use crate::coord::ellipsoid::{Ellipsoid, EllipsoidType};

const MAX_TERMS : usize = 8;


#[derive(Copy, Clone)]
pub struct TransMercCoeffs {
    pub n1       : f64,
    pub aCoeff   : [f64; 8],
    pub bCoeff   : [f64; 8],
    pub R4oa     : f64,
}

impl TransMercCoeffs {
    pub fn new() -> TransMercCoeffs {
        TransMercCoeffs {
            n1: 0.0,
            aCoeff: [0.0; MAX_TERMS],
            bCoeff: [0.0; MAX_TERMS],
            R4oa: 0.0,
        }
    }

    pub fn get_a_coeffs( ellipsoid: Ellipsoid ) -> [f64; MAX_TERMS] {

        let mut res: [f64; MAX_TERMS] = [0.0; MAX_TERMS];

        // Misc coeffients used in computations
        let mut coeff: f64 = 0.0;
        let n1: f64 = 1.0 / (2.0 * ellipsoid.inv_f - 1.0);
        let n2: f64 = n1 * n1;
        let n3: f64 = n2 * n1;
        let n4: f64 = n3 * n1;
        let n5: f64 = n4 * n1;
        let n6: f64 = n5 * n1;
        let n7: f64 = n6 * n1;
        let n8: f64 = n7 * n1;
        let n9: f64 = n8 * n1;
        let n10: f64 = n9 * n1;

        match ellipsoid.code {

            // World
            EllipsoidType::WGS84 => {
                res[0] = 8.3773182062446983032e-04;
                res[1] = 7.608527773572489156e-07;
                res[2] = 1.19764550324249210e-09;
                res[3] = 2.4291706803973131e-12;
                res[4] = 5.711818369154105e-15;
                res[5] = 1.47999802705262e-17;
                res[6] = 0.0; // LSC-16618
                res[7] = 0.0; // LSC-16618
                return res;
            },

            // Some Custom Ellipsoid
            EllipsoidType::CUSTOM => {

                // computation below is for user defined ellipsoids
                // Computation of coefficient a2
                coeff = 0.0;
                coeff += (-18975107.0) * n8 / 50803200.0;
                coeff += (72161.0) * n7 / 387072.0;
                coeff += (7891.0) * n6 / 37800.0;
                coeff += (-127.0) * n5 / 288.0;
                coeff += (41.0) * n4 / 180.0;
                coeff += (5.0) * n3 / 16.0;
                coeff += (-2.0) * n2 / 3.0;
                coeff += (1.0) * n1 / 2.0;
                res[0] = coeff;

                //   Computation of coefficient a4
                coeff = 0.0;
                coeff += (148003883.0) * n8 / 174182400.0;
                coeff += (13769.0) * n7 / 28800.0;
                coeff += (-1983433.0) * n6 / 1935360.0;
                coeff += (281.0) * n5 / 630.0;
                coeff += (557.0) * n4 / 1440.0;
                coeff += (-3.0) * n3 / 5.0;
                coeff += (13.0) * n2 / 48.0;
                res[1] = coeff;

                //   Computation of coefficient a6
                coeff = 0.0;
                coeff += (79682431.0) * n8 / 79833600.0;
                coeff += (-67102379.0) * n7 / 29030400.0;
                coeff += (167603.0) * n6 / 181440.0;
                coeff += (15061.0) * n5 / 26880.0;
                coeff += (-103.0) * n4 / 140.0;
                coeff += (61.0) * n3 / 240.0;
                res[2] = coeff;

                //   Computation of coefficient a8
                coeff = 0.0;
                coeff += (-40176129013.0) * n8 / 7664025600.0;
                coeff += (97445.0) * n7 / 49896.0;
                coeff += (6601661.0) * n6 / 7257600.0;
                coeff += (-179.0) * n5 / 168.0;
                coeff += (49561.0) * n4 / 161280.0;
                res[3] = coeff;

                //   Computation of coefficient a10
                coeff = 0.0;
                coeff += (2605413599.0) * n8 / 622702080.0;
                coeff += (14644087.0) * n7 / 9123840.0;
                coeff += (-3418889.0) * n6 / 1995840.0;
                coeff += (34729.0) * n5 / 80640.0;
                res[4] = coeff;

                //   Computation of coefficient a12
                coeff = 0.0;
                coeff += (175214326799.0) * n8 / 58118860800.0;
                coeff += (-30705481.0) * n7 / 10378368.0;
                coeff += (212378941.0) * n6 / 319334400.0;
                res[5] = coeff;

                //   Computation of coefficient a14
                coeff = 0.0;
                coeff += (-16759934899.0) * n8 / 3113510400.0;
                coeff += (1522256789.0) * n7 / 1383782400.0;
                res[6] = coeff;

                //   Computation of coefficient a16
                coeff = 0.0;
                coeff += (1424729850961.0) * n8 / 743921418240.0;
                res[7] = coeff;
                return res;
            }
        }
    }

    pub fn get_b_coeffs( ellipsoid: Ellipsoid ) -> [f64; MAX_TERMS] {

        let mut res: [f64; MAX_TERMS] = [0.0; MAX_TERMS];

        // Misc coeffients used in computations
        let mut coeff: f64 = 0.0;
        let n1: f64 = 1.0 / (2.0 * ellipsoid.inv_f - 1.0);
        let n2: f64 = n1 * n1;
        let n3: f64 = n2 * n1;
        let n4: f64 = n3 * n1;
        let n5: f64 = n4 * n1;
        let n6: f64 = n5 * n1;
        let n7: f64 = n6 * n1;
        let n8: f64 = n7 * n1;
        let n9: f64 = n8 * n1;
        let n10: f64 = n9 * n1;

        match ellipsoid.code {
            EllipsoidType::WGS84 => {
                res[0] = -8.3773216405794867707e-04;
                res[1] = -5.905870152220365181e-08;
                res[2] = -1.67348266534382493e-10;
                res[3] = -2.1647981104903862e-13;
                res[4] = -3.787930968839601e-16;
                res[5] = -7.23676928796690e-19;
                res[6] = 0.0; // LSC-16618
                res[7] = 0.0; // LSC-16618
                return res;
            },

            EllipsoidType::CUSTOM => {
                //   Computation of coefficient b2
                coeff = 0.0;
                coeff += (-7944359.0) * n8 / 67737600.0;
                coeff += (5406467.0) * n7 / 38707200.0;
                coeff += (-96199.0) * n6 / 604800.0;
                coeff += (81.0) * n5 / 512.0;
                coeff += (1.0) * n4 / 360.0;
                coeff += (-37.0) * n3 / 96.0;
                coeff += (2.0) * n2 / 3.0;
                coeff += (-1.0) * n1 / 2.0;
                res[0] = coeff;

                //   Computation of coefficient b4
                coeff = 0.0;
                coeff += (-24749483.0) * n8 / 348364800.0;
                coeff += (-51841.0) * n7 / 1209600.0;
                coeff += (1118711.0) * n6 / 3870720.0;
                coeff += (-46.0) * n5 / 105.0;
                coeff += (437.0) * n4 / 1440.0;
                coeff += (-1.0) * n3 / 15.0;
                coeff += (-1.0) * n2 / 48.0;
                res[1] = coeff;

                //   Computation of coefficient b6
                coeff = 0.0;
                coeff += (6457463.0) * n8 / 17740800.0;
                coeff += (-9261899.0) * n7 / 58060800.0;
                coeff += (-5569.0) * n6 / 90720.0;
                coeff += (209.0) * n5 / 4480.0;
                coeff += (37.0) * n4 / 840.0;
                coeff += (-17.0) * n3 / 480.0;
                res[2] = coeff;

                //   Computation of coefficient b8
                coeff = 0.0;
                coeff += (-324154477.0) * n8 / 7664025600.0;
                coeff += (-466511.0) * n7 / 2494800.0;
                coeff += (830251.0) * n6 / 7257600.0;
                coeff += (11.0) * n5 / 504.0;
                coeff += (-4397.0) * n4 / 161280.0;
                res[3] = coeff;

                //   Computation of coefficient b10
                coeff = 0.0;
                coeff += (-22894433.0) * n8 / 124540416.0;
                coeff += (8005831.0) * n7 / 63866880.0;
                coeff += (108847.0) * n6 / 3991680.0;
                coeff += (-4583.0) * n5 / 161280.0;
                res[4] = coeff;

                //   Computation of coefficient b12
                coeff = 0.0;
                coeff += (2204645983.0) * n8 / 12915302400.0;
                coeff += (16363163.0) * n7 / 518918400.0;
                coeff += (-20648693.0) * n6 / 638668800.0;
                res[5] = coeff;

                //   Computation of coefficient b14
                coeff = 0.0;
                coeff += (497323811.0) * n8 / 12454041600.0;
                coeff += (-219941297.0) * n7 / 5535129600.0;
                res[6] = coeff;

                //   Computation of coefficient b16
                coeff = 0.0;
                coeff += (-191773887257.0) * n8 / 3719607091200.0;
                res[7] = coeff;
                return res;
            },
        }
    }

    // This method is derived from GeoTrans, file TransverseMercator.cpp
    // Algorithm developed by: C. Rollins   April 18, 2006
    // INPUT
    // -----
    // inv_fla    Inverse flattening (reciprocal flattening)
    // ellipsoid  Ellipsoidal parameters
    //
    // OUTPUT
    //     ------
    //        n1        Helmert's "n"
    //        aCoeff    Coefficients for omega as a trig series in chi
    //        bBoeff    Coefficients for chi as a trig series in omega
    //        R4oa      Ratio "R4 over a", i.e. R4/a
    //
    //     EXPLANATIONS
    //     ------------
    //        omega is rectifying latitude
    //        chi is conformal latitude
    //        psi is geocentric latitude
    //        phi is geodetic latitude, commonly, "the latitude"
    //        R4 is the meridional isoperimetric radius
    //        "a" is the semi-major axis of the ellipsoid
    //        "b" is the semi-minor axis of the ellipsoid
    //        Helmert's n = (a - b)/(a + b)
    //
    //        This calculation depends only on the shape of the ellipsoid and are
    //        independent of the ellipsoid size.
    //
    //        The array Acoeff(8) stores eight coefficients corresponding
    //           to k = 2, 4, 6, 8, 10, 12, 14, 16 in the notation "a sub k".
    //        Likewise Bcoeff(8) etc.
    pub fn generate_coefficents( ellipsoid: Ellipsoid ) -> TransMercCoeffs {

        // checks ellipsoid code and assigns values for corresponding coefficients.
        // Uses default computation if ellipsoid code isn't found. This will be
        // for user defined ellipsoids.

        //  Generate Coefficients for Transverse Mercator algorithms
        // double n2, n3, n4, n5, n6, n7, n8, n9, n10, coeff;

        let mut res = TransMercCoeffs::new();

        // Set the Coefficients
        res.aCoeff = TransMercCoeffs::get_a_coeffs( ellipsoid );
        res.bCoeff = TransMercCoeffs::get_b_coeffs( ellipsoid );

        //  Create base coefficients
        let mut coeff: f64 = 0.0;
        let n1: f64 = 1.0 / (2.0 * ellipsoid.inv_f - 1.0);
        let n2: f64 = n1 * n1;
        let n3: f64 = n2 * n1;
        let n4: f64 = n3 * n1;
        let n5: f64 = n4 * n1;
        let n6: f64 = n5 * n1;
        let n7: f64 = n6 * n1;
        let n8: f64 = n7 * n1;
        let n9: f64 = n8 * n1;
        let n10: f64 = n9 * n1;

        coeff += 49.0 * n10 / 65536.0;
        coeff += 25.0 * n8 / 16384.0;
        coeff +=      n6 / 256.0;
        coeff +=      n4 / 64.0;
        coeff +=      n2 / 4.0;
        coeff += 1.0;
        res.R4oa = coeff / (1.0 + n1);

        return res;
    }
}

// if (( strcmp( ellipsoidCode, "AA") == 0) || (strcmp(ellipsoidCode, "AM") == 0))
// {
// aCoeff[0] = 8.3474517669594013740e-04;
// aCoeff[1] = 7.554352936725572895e-07;
// aCoeff[2] = 1.18487391005135489e-09;
// aCoeff[3] = 2.3946872955703565e-12;
// aCoeff[4] = 5.610633978440270e-15;
// aCoeff[5] = 1.44858956458553e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3474551646761162264e-04;
// bCoeff[1] = -5.863630361809676570e-08;
// bCoeff[2] = -1.65562038746920803e-10;
// bCoeff[3] = -2.1340335537652749e-13;
// bCoeff[4] = -3.720760760132477e-16;
// bCoeff[5] = -7.08304328877781e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if (( strcmp( ellipsoidCode, "EA") == 0) || ( strcmp( ellipsoidCode, "EB") == 0) ||
// ( strcmp( ellipsoidCode, "EC") == 0) || ( strcmp( ellipsoidCode, "ED") == 0) ||
// ( strcmp( ellipsoidCode, "EE") == 0))
// aCoeff[0] = 8.3064943111192510534e-04;
// aCoeff[1] = 7.480375027595025021e-07;
// aCoeff[2] = 1.16750772278215999e-09;
// aCoeff[3] = 2.3479972304395461e-12;
// aCoeff[4] = 5.474212231879573e-15;
// aCoeff[5] = 1.40642257446745e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3064976590443772201e-04;
// bCoeff[1] = -5.805953517555717859e-08;
// bCoeff[2] = -1.63133251663416522e-10;
// bCoeff[3] = -2.0923797199593389e-13;
// bCoeff[4] = -3.630200927775259e-16;
// bCoeff[5] = -6.87666654919219e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if (( strcmp( ellipsoidCode, "BN") == 0) || ( strcmp( ellipsoidCode, "BR") == 0))
// aCoeff[0] = 8.3522527226849818552e-04;
// aCoeff[1] = 7.563048340614894422e-07;
// aCoeff[2] = 1.18692075307408346e-09;
// aCoeff[3] = 2.4002054791393298e-12;
// aCoeff[4] = 5.626801597980756e-15;
// aCoeff[5] = 1.45360057224474e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3522561262703079182e-04;
// bCoeff[1] = -5.870409978661008580e-08;
// bCoeff[2] = -1.65848307463131468e-10;
// bCoeff[3] = -2.1389565927064571e-13;
// bCoeff[4] = -3.731493368666479e-16;
// bCoeff[5] = -7.10756898071999e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if (( strcmp( ellipsoidCode, "KA") == 0) || ( strcmp( ellipsoidCode, "HE") == 0) ||
// ( strcmp( ellipsoidCode, "FA") == 0))
// aCoeff[0] = 8.3761175713442343106e-04;
// aCoeff[1] = 7.606346200814720197e-07;
// aCoeff[2] = 1.19713032035541037e-09;
// aCoeff[3] = 2.4277772986483520e-12;
// aCoeff[4] = 5.707722772225013e-15;
// aCoeff[5] = 1.47872454335773e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3761210042019176501e-04;
// bCoeff[1] = -5.904169154078546237e-08;
// bCoeff[2] = -1.67276212891429215e-10;
// bCoeff[3] = -2.1635549847939549e-13;
// bCoeff[4] = -3.785212121016612e-16;
// bCoeff[5] = -7.23053625983667e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if ( strcmp( ellipsoidCode, "WD") == 0)
// aCoeff[0] = 8.3772481044362217923e-04;
// aCoeff[1] = 7.608400388863560936e-07;
// aCoeff[2] = 1.19761541904924067e-09;
// aCoeff[3] = 2.4290893081322466e-12;
// aCoeff[4] = 5.711579173743133e-15;
// aCoeff[5] = 1.47992364667635e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3772515386847544554e-04;
// bCoeff[1] = -5.905770828762463028e-08;
// bCoeff[2] = -1.67344058948464124e-10;
// bCoeff[3] = -2.1647255130188214e-13;
// bCoeff[4] = -3.787772179729998e-16;
// bCoeff[5] = -7.23640523525528e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if ( strcmp( ellipsoidCode, "RF") == 0)
// aCoeff[0] = 8.3773182472855134012e-04;
// aCoeff[1] = 7.608527848149655006e-07;
// aCoeff[2] = 1.19764552085530681e-09;
// aCoeff[3] = 2.4291707280369697e-12;
// aCoeff[4] = 5.711818509192422e-15;
// aCoeff[5] = 1.47999807059922e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3773216816203523672e-04;
// bCoeff[1] = -5.905870210369121594e-08;
// bCoeff[2] = -1.67348268997717031e-10;
// bCoeff[3] = -2.1647981529928124e-13;
// bCoeff[4] = -3.787931061803592e-16;
// bCoeff[5] = -7.23676950110361e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if (( strcmp( ellipsoidCode, "SA") == 0) || ( strcmp( ellipsoidCode, "AN") == 0))
// aCoeff[0] = 8.3775209887947194075e-04;
// aCoeff[1] = 7.608896263599627157e-07;
// aCoeff[2] = 1.19773253021831769e-09;
// aCoeff[3] = 2.4294060763606098e-12;
// aCoeff[4] = 5.712510331613028e-15;
// aCoeff[5] = 1.48021320370432e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3775244233790270051e-04;
// bCoeff[1] = -5.906157468586898015e-08;
// bCoeff[2] = -1.67360438158764851e-10;
// bCoeff[3] = -2.1650081225048788e-13;
// bCoeff[4] = -3.788390325953455e-16;
// bCoeff[5] = -7.23782246429908e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

//  else if ( strcmp( ellipsoidCode, "ID") == 0)
// aCoeff[0] = 8.3776052087969078729e-04;
// aCoeff[1] = 7.609049308144604484e-07;
// aCoeff[2] = 1.19776867565343872e-09;
// aCoeff[3] = 2.4295038464530901e-12;
// aCoeff[4] = 5.712797738386076e-15;
// aCoeff[5] = 1.48030257891140e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.3776086434848497443e-04;
// bCoeff[1] = -5.906276799395007586e-08;
// bCoeff[2] = -1.67365493472742884e-10;
// bCoeff[3] = -2.1650953495573773e-13;
// bCoeff[4] = -3.788581120060625e-16;
// bCoeff[5] = -7.23825990889693e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

//  else if (( strcmp( ellipsoidCode, "IN") == 0) || ( strcmp( ellipsoidCode, "HO") == 0))
// aCoeff[0] = 8.4127599100356448089e-04;
// aCoeff[1] = 7.673066923431950296e-07;
// aCoeff[2] = 1.21291995794281190e-09;
// aCoeff[3] = 2.4705731165688123e-12;
// aCoeff[4] = 5.833780550286833e-15;
// aCoeff[5] = 1.51800420867708e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618
// bCoeff[0] = -8.4127633881644851945e-04;
// bCoeff[1] = -5.956193574768780571e-08;
// bCoeff[2] = -1.69484573979154433e-10;
// bCoeff[3] = -2.2017363465021880e-13;
// bCoeff[4] = -3.868896221495780e-16;
// bCoeff[5] = -7.42279219864412e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if ( strcmp( ellipsoidCode, "WO") == 0)
// aCoeff[0] = 8.4411652150600103279e-04;
// aCoeff[1] = 7.724989750172583427e-07;
// aCoeff[2] = 1.22525529789972041e-09;
// aCoeff[3] = 2.5041361775549209e-12;
// aCoeff[4] = 5.933026083631383e-15;
// aCoeff[5] = 1.54904908794521e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.4411687285559594196e-04;
// bCoeff[1] = -5.996681687064322548e-08;
// bCoeff[2] = -1.71209836918814857e-10;
// bCoeff[3] = -2.2316811233502163e-13;
// bCoeff[4] = -3.934782433323038e-16;
// bCoeff[5] = -7.57474665717687e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if ( strcmp( ellipsoidCode, "CC") == 0)
// aCoeff[0] = 8.4703742793654652315e-04;
// aCoeff[1] = 7.778564517658115212e-07;
// aCoeff[2] = 1.23802665917879731e-09;
// aCoeff[3] = 2.5390045684252928e-12;
// aCoeff[4] = 6.036484469753319e-15;
// aCoeff[5] = 1.58152259295850e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.4703778294785813001e-04;
// bCoeff[1] = -6.038459874600183555e-08;
// bCoeff[2] = -1.72996106059227725e-10;
// bCoeff[3] = -2.2627911073545072e-13;
// bCoeff[4] = -4.003466873888566e-16;
// bCoeff[5] = -7.73369749524777e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if ( strcmp( ellipsoidCode, "CG") == 0)
// aCoeff[0] = 8.5140099460764136776e-04;
// aCoeff[1] = 7.858945456038187774e-07;
// aCoeff[2] = 1.25727085106103462e-09;
// aCoeff[3] = 2.5917718627340128e-12;
// aCoeff[4] = 6.193726879043722e-15;
// aCoeff[5] = 1.63109098395549e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618
// bCoeff[0] = -8.5140135513650084564e-04;
// bCoeff[1] = -6.101145475063033499e-08;
// bCoeff[2] = -1.75687742410879760e-10;
// bCoeff[3] = -2.3098718484594067e-13;
// bCoeff[4] = -4.107860472919190e-16;
// bCoeff[5] = -7.97633133452512e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618

// else if ( strcmp( ellipsoidCode, "CD") == 0)
// aCoeff[0] = 8.5140395445291970541e-04;
// aCoeff[1] = 7.859000119464140978e-07;
// aCoeff[2] = 1.25728397182445579e-09;
// aCoeff[3] = 2.5918079321459932e-12;
// aCoeff[4] = 6.193834639108787e-15;
// aCoeff[5] = 1.63112504092335e-17;
// aCoeff[6] = 0.0; // LSC-16618
// aCoeff[7] = 0.0; // LSC-16618

// bCoeff[0] = -8.5140431498554106268e-04;
// bCoeff[1] = -6.101188106187092184e-08;
// bCoeff[2] = -1.75689577596504470e-10;
// bCoeff[3] = -2.3099040312610703e-13;
// bCoeff[4] = -4.107932016207395e-16;
// bCoeff[5] = -7.97649804397335e-19;
// bCoeff[6] = 0.0; // LSC-16618
// bCoeff[7] = 0.0; // LSC-16618