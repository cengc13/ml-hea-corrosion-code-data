from qmpy.analysis.miedema import Miedema

comp = {'Fe':1, 'Se':1}

obj = Miedema(comp)

print(obj.H_mix)