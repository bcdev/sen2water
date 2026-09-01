# -*- coding: utf-8 -*-

"""..."""

__author__ = "Martin Böttcher, Olaf Danne (Brockmann Consult GmbH)"
__copyright__ = "Copyright 2023-2026, Brockmann Consult GmbH"
__license__ = "MIT"
__version__ = "0.6"
__email__ = "info@brockmann-consult.de"
__status__ = "Development"

import os

# changes in 1.1:
# ...

import numpy as np
import pytest

import tensorflow as tf


# Test centre pixel px = 1300, py = 1200 from S2A_MSIL1C_20170815T102021_N0205_R065_T33UUA_20170815T102513.SAFE,
# resampled to 60m. Consider also 3x3 box around.

# Definition of input in SNAP C2rccMsiAlgorithm.java:

#  //  (9.2) compute angles
# double cos_sun = cos(toRadians(sun_zeni));
# double cos_view = cos(toRadians(view_zeni));
# double sin_view = sin(toRadians(view_zeni));
#
# double cos_azi_diff = cos(toRadians(view_azi - sun_azi));
# double azi_diff_rad = acos(cos_azi_diff);
# double sin_azi_diff = sin(azi_diff_rad);
# double azi_diff_deg = toDegrees(azi_diff_rad);
#
# double x = sin_view * cos_azi_diff;
# double y = sin_view * sin_azi_diff;
# double z = cos_view;

# // (9.4) set input to all atmosphere NNs
# //nn_in=[sun_zeni,x,y,z,temperature, salinity, alti_press, log_rtosa];
# int ancNnInputCount = 7;
# double[] nn_in = new double[ancNnInputCount + log_rtosa.length];
# nn_in[0] = sun_zeni;
# nn_in[1] = x;  # = sin_view * cos_azi_diff =
# nn_in[2] = y;
# nn_in[3] = z;  # = cos_view = cos(toRadians(view_zeni));
# nn_in[4] = temperature;
# nn_in[5] = salinity;
# nn_in[6] = alti_press;
# System.arraycopy(log_rtosa, 0, nn_in, ancNnInputCount, log_rtosa.length);

nn_in_1299_1200 = np.array(
    [41.49864196777344, -0.01205048826557415, 0.01869579483574073, 0.9997525958896136, 15.0, 35.0, 1000.0,
     -2.0537660112576814, -2.305917292115702, -2.4931280685278905, -2.9063789330443397, -3.035973302453955,
     -3.2206552364722927, -3.2538968664513894, -3.4067620091879363], dtype=np.float32)

nn_in_1300_1200 = np.array([41.49850082397461, -0.012114469857435216, 0.018628545618525064, 0.9997530779688611,
                            15.0, 35.0, 1000.0, -2.0475306092209555, -2.2947666674186156, -2.4839837727443017,
                            -2.8875460408637563, -3.0212826092254876, -3.1957163001382436, -3.225653662563658,
                            -3.356677906911101], dtype=np.float32)
nn_in_1301_1200 = np.array(
    [41.49835968017578, -0.012178101939365298, 0.01856108185165093, 0.9997535596704074, 15.0, 35.0, 1000.0,
     -2.0467539141328728, -2.2897390277613643, -2.4749223400882476, -2.867231455658479, -2.9844686582816973,
     -3.1737909515210676, -3.173854175628189, -3.3117265211691613], dtype=np.float32)

nn_in_1299_1199 = np.array(
    [41.499168395996094, -0.012035523322517012, 0.018718930792967072, 0.9997523432372245, 15.0, 35.0, 1000.0,
     -2.0592540559041974, -2.312052205308307, -2.5103352542519675, -2.927513061276227, -3.061670575327375,
     -3.261898155675378, -3.2829609718098625, -3.4658771898127965], dtype=np.float32)

nn_in_1300_1199 = np.array(
    [41.499027252197266, -0.012099478436432541, 0.018651696344171798, 0.999752827875496, 15.0, 35.0, 1000.0,
     -2.056114339883183, -2.3008333447438463, -2.4957559226160164, -2.9082817555445186, -3.042335945818706,
     -3.2131080035284914, -3.25129606559811, -3.394749847919438], dtype=np.float32)
nn_in_1301_1199 = np.array(
    [41.49888610839844, -0.012163084830578897, 0.018584259773197942, 0.999753311900534, 15.0, 35.0, 1000.0,
     -2.0529844508571644, -2.3008333677798665, -2.4957560407211337, -2.8988030677756456, -3.019201392611685,
     -3.215617427154215, -3.235831848331168, -3.379934764256126], dtype=np.float32)

nn_in_1299_1201 = np.array(
    [41.49811553955078, -0.012065413610398512, 0.01867265565793312, 0.9997528483205688, 15.0, 35.0, 1000.0,
     -2.047530618272058, -2.299819787619027, -2.494441652780426, -2.8988033230695125, -3.0233681661534972,
     -3.2106049399479466, -3.2332776068328233, -3.3682387341340014], dtype=np.float32)

nn_in_1300_1201 = np.array(
    [41.49797439575195, -0.012129428279344332, 0.018605377351934998, 0.9997533280280733, 15.0, 35.0, 1000.0,
     -2.046753923181568, -2.293759186013127, -2.4865883272249882, -2.887546224656819, -3.0109198947924396,
     -3.190802324656291, -3.210578027736293, -3.3595556136795848], dtype=np.float32)
nn_in_1301_1201 = np.array(
    [41.497833251953125, -0.012193087176950798, 0.018537885478991947, 0.9997538074081354, 15.0, 35.0, 1000.0,
     -2.045977830883952, -2.295775338731734, -2.4826847316825083, -2.8782611629873047, -3.0027063184182805,
     -3.1762035372384263, -3.1908243074187705, -3.3200024374123895], dtype=np.float32)

nn_expected_result_rtosa_trans_1299_1200 = np.array(
    [0.8413484524781636, 0.8777395802934942, 0.9159607714939266, 0.9431746274141641, 0.9484070118147985,
      0.9527421907744882, 0.9567231234765484, 0.9605067045532839, 0.8725681711429314, 0.9039329673592331,
      0.9360685917958538, 0.95840325799782, 0.9626217533863966, 0.9660593914622994, 0.9691799860038679,
      0.972142906035206])
nn_expected_result_rtosa_trans_1300_1200 = np.array(
    [0.8408185626136915, 0.8772525906338025, 0.915509626091523, 0.9427289328411357, 0.9479384825989265,
      0.9522327900235773, 0.9561653423992285, 0.9598742761401144, 0.8724045450248956, 0.9037699831350547,
      0.9358993530953172, 0.9582073197097325, 0.9623984260329231, 0.9658014841207325, 0.9688743646612192,
      0.9717736374959793])
nn_expected_result_rtosa_trans_1301_1200 = np.array(
    [0.8424705881585794, 0.8788682720476235, 0.9169959729931921, 0.9440141017708938, 0.9491596682347101,
      0.9534012813763244, 0.9572837139415608, 0.9609423329375477, 0.8734678459247889, 0.9047869770397077,
      0.9368293307587778, 0.9590266874335949, 0.9631921299028506, 0.9665742569169816, 0.9696369535651419,
      0.9725219916222886])

nn_expected_result_rtosa_trans_1299_1199 = np.array(
    [0.8434842553010391, 0.879932570082928, 0.918151713170703, 0.9452962386763677, 0.950558175855672,
      0.9549574816543827, 0.9590204143613271, 0.9629488576977661, 0.8735257714506244, 0.9049462459171514,
      0.9371513889352723, 0.9595704119108576, 0.9638359291520906, 0.9673376264832892, 0.9705298312067363,
      0.9736683648670035])
nn_expected_result_rtosa_trans_1300_1199 = np.array(
    [0.8430411360667944, 0.8795343948309978, 0.9177793910667084, 0.9448734711750256, 0.9501036523736499,
      0.9544443220156149, 0.9584439160769421, 0.9622913985468368, 0.8734500090177534, 0.9048627040197618,
      0.9370327183752851, 0.9594035289184774, 0.9636396763103003, 0.9671083073351044, 0.9702627488122285,
      0.9733321998109775])
nn_expected_result_rtosa_trans_1301_1199 = np.array(
    [0.8421680471796189, 0.8786109408809799, 0.9168363018707559, 0.9439826163049978, 0.9491837058841586,
      0.9534879035691317, 0.9574446774667325, 0.9611989826297069, 0.8734277038444954, 0.9047630353697129,
      0.936832847571884, 0.9590821479995758, 0.9632726874189027, 0.9666842073123266, 0.969782155032275,
      0.9727267231353334])

nn_expected_result_rtosa_trans_1299_1201 = np.array(
    [0.8400063611680165, 0.8765059451506938, 0.914873098240069, 0.9422325055214493, 0.9474706847005192,
      0.9517769633316724, 0.9557327999468285, 0.9594655431642831, 0.8723908864953142, 0.9037495691085461,
      0.9358668573684762, 0.9581501980098679, 0.9623307927080569, 0.9657150174879742, 0.9687766238728803,
      0.971629981245365])
nn_expected_result_rtosa_trans_1300_1201 = np.array(
    [0.8415594211917907, 0.8779352459188621, 0.9161231667681083, 0.9432644795556717, 0.9484804862980586,
      0.952800984781812, 0.9567570022878111, 0.9605067045532839, 0.872522770648488, 0.9038882048832757,
      0.9360254555759326, 0.9583622788192409, 0.9625760013819138, 0.9660091720660644, 0.969123089002441,
      0.972085048035166])
nn_expected_result_rtosa_trans_1301_1201 = np.array(
    [0.8403814171341635, 0.8767904488431573, 0.9150546109963271, 0.9423047164413895, 0.9475020399905825,
      0.9517709300277176, 0.9556747540191808, 0.9593397505477533, 0.8723999925747089, 0.9037373163257666,
      0.9358307177557059, 0.9580802383243914, 0.9622444937935538, 0.965613326857974, 0.9686554695833226,
      0.9714768135594638])


tf_model_root_path = (os.getcwd() + os.sep + "sen2water" + os.sep + "msic2rcc" + os.sep + "resources" + os.sep +
                      "c2rcc-msi-nets-tensorflow" + os.sep + "msi" + os.sep + "complex_20200321")


def test_nn_rtosa_trans():
    # C2rccMsiAlgorithm line 373
    tf_model_path = tf_model_root_path + os.sep + "rtosa_trans" + os.sep + "55x55x55_107567.4.tf"

    model = tf.saved_model.load(tf_model_path)
    infer = model.signatures["serving_default"]

    nn_result = infer(inputs=tf.convert_to_tensor(nn_in_1300_1200))["predictions"].numpy()
    # print("nn_result: " + str(nn_result))
    np.testing.assert_allclose(nn_result[0], nn_expected_result_rtosa_trans_1300_1200, rtol=1.E-5, atol=1.E-5)

    nn_result = infer(inputs=tf.convert_to_tensor(nn_in_1299_1200))["predictions"].numpy()
    # print("nn_result: " + str(nn_result))
    np.testing.assert_allclose(nn_result[0], nn_expected_result_rtosa_trans_1299_1200, rtol=1.E-5, atol=1.E-5)

    nn_result = infer(inputs=tf.convert_to_tensor(nn_in_1301_1200))["predictions"].numpy()
    # print("nn_result: " + str(nn_result))
    np.testing.assert_allclose(nn_result[0], nn_expected_result_rtosa_trans_1301_1200, rtol=1.E-5, atol=1.E-5)

    print("done test_nn_rtosa_trans")


def test_nn_rtosa_trans_array_3x1():
    # C2rccMsiAlgorithm line 373
    tf_model_path = tf_model_root_path + os.sep + "rtosa_trans" + os.sep + "55x55x55_107567.4.tf"

    model = tf.saved_model.load(tf_model_path)
    infer = model.signatures["serving_default"]

    nn_input_array = np.array([nn_in_1299_1200, nn_in_1300_1200, nn_in_1301_1200])
    nn_result = infer(inputs=tf.convert_to_tensor(nn_input_array))["predictions"].numpy()
    print("nn_result: " + str(nn_result))
    nn_expected_result_array = np.array(
        [nn_expected_result_rtosa_trans_1299_1200, nn_expected_result_rtosa_trans_1300_1200,
         nn_expected_result_rtosa_trans_1301_1200])
    np.testing.assert_allclose(nn_result, nn_expected_result_array, rtol=1.E-5, atol=1.E-5)

    print("done test_nn_rtosa_trans_array_3x1")

def test_nn_rtosa_trans_array_3x3():
    # C2rccMsiAlgorithm line 373
    tf_model_path = tf_model_root_path + os.sep + "rtosa_trans" + os.sep + "55x55x55_107567.4.tf"

    model = tf.saved_model.load(tf_model_path)
    infer = model.signatures["serving_default"]

    nn_input_array = np.array([nn_in_1299_1199, nn_in_1300_1199, nn_in_1301_1199, nn_in_1299_1200, nn_in_1300_1200, nn_in_1301_1200, nn_in_1299_1201, nn_in_1300_1201, nn_in_1301_1201])
    nn_result = infer(inputs=tf.convert_to_tensor(nn_input_array))["predictions"].numpy()
    print("nn_result: " + str(nn_result))
    nn_expected_result_array = np.array(
        [nn_expected_result_rtosa_trans_1299_1199, nn_expected_result_rtosa_trans_1300_1199, nn_expected_result_rtosa_trans_1301_1199, nn_expected_result_rtosa_trans_1299_1200, nn_expected_result_rtosa_trans_1300_1200, nn_expected_result_rtosa_trans_1301_1200, nn_expected_result_rtosa_trans_1299_1201, nn_expected_result_rtosa_trans_1300_1201, nn_expected_result_rtosa_trans_1301_1201])
    np.testing.assert_allclose(nn_result, nn_expected_result_array, rtol=1.E-5, atol=1.E-5)

    print("done test_nn_rtosa_trans_array_3x3")


@pytest.mark.skip
def test_nn_iop_rw_random():
    # TODO

    tf_model_path = tf_model_root_path + os.sep + "iop_rw" + os.sep + "55x55x55_2775.5.tf"

    model = tf.saved_model.load(tf_model_path)
    infer = model.signatures["serving_default"]

    # n_features = 12  # use 13 if this model includes B10
    n_features = 10  # use 13 if this model includes B10

    x = np.random.rand(1, n_features).astype(np.float32)
    y = infer(inputs=tf.convert_to_tensor(x))["predictions"].numpy()

    print("y: " + str(y))

    print("done test_nn_iop_rw_random")
