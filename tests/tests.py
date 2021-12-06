print("Testing imports ....")

try:
	print("Imporintg simunet.analysis.network")
	from simunet.analysis.network import generate_subnetworks
	print("Importing simunet.analysis.methods")
	from simunet.analysis.methods import *
	print("Importing simunet.analysis.scores")
	from simunet.analysis.scores import *
	print("Importing simunet.common.errors")
	from simunet.common.errors import *
	print("Importing simunet.common.parsers")
	from simunet.common.parsers import *
	print("Importing simunet.common.common")
	from simunet.common.common import *
except ImportError:
    print("ERROR: last import Failed")
finally:
	print("Testing complete")