#import shablona as sb
import LHconnectivity as lh


def test_average_structurals():
    WD='/home/share/LeftHand/LHconnectivity/Stockholm'
    subject = 'LH1001'

    lh.average_structurals(WD, subject)


test_average_structurals()
