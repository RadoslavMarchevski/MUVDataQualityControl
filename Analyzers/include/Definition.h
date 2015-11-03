//
//  Definition.h
//
//
//  Created by Riccardo Aliberti on 03/05/13.
//
//
#ifndef _Value_
#define _Value_

//PDG Values
const double MuMass   = 0.10565837;// [ GeV ]
const double Pi0Mass  = 0.1349766; // [ GeV ]
const double PiPlMass = 0.13957;   // [ GeV ]
const double KMass = 0.493677;     // [ GeV ]

//Beam Parameters
const double KEnergy = 74.8; // [ GeV ]
const double XAngle = 0.0012; // [ rad ]

//Detector Parameters
const double Ztrim = 101.8; // [ M ]
const double ZMagnetStart = 196.345; // [ M ]
const double ZMagnetEnd = 197.645; // [ M ]
const double ZMagnetMean = (ZMagnetEnd + ZMagnetStart)/2.;
const double MagnetKick = 0.2574; // [ M ]
//const double MagnetKick = 0.300; // [ M ]

const double ZCedarStart = 69.165; // [ M ]
const double ZCedarEnd = 79.440; // [ M ]

const double ZRICHStart= 219.385; // [ M ]
const double ZCHODStart = 238.960; // [ M ]
const double RCHODmin = 0.15; // [ M ]
const double RCHODmax = 1.200; // [ M ]

const double ZLAV1_Start = 121.953; // [ M ]
const double ZLAV2_Start = 129.563; // [ M ]
const double ZLAV3_Start = 137.173; // [ M ]
const double ZLAV4_Start = 144.783; // [ M ]
const double ZLAV5_Start = 152.393; // [ M ]
const double ZLAV6_Start = 165.903; // [ M ]
const double ZLAV7_Start = 173.413; // [ M ]
const double ZLAV8_Start = 180.923; // [ M ]
const double ZLAV9_Start = 193.184; // [ M ]
const double ZLAV10_Start = 203.642;// [ M ]
const double ZLAV11_Start = 218.203;// [ M ]
const double ZLAV12_Start = 238.835;// [ M ]


const double ZLKrStart  = 241.495; // [ M ]
const double ZMUV1Start = 244.341;
const double XMUV1max = 1.200;
const double YMUV1max = 1.200;
const double RMUV1max = 1.200;
const double RMUV1min = 0.15;

const double ZMUV2Start = 245.290;
const double XMUV2max = 1.200;
const double YMUV2max = 1.200;
const double RMUV2max = 1.200;
const double RMUV2min = 0.15;


const double ZMUV3Start = 246.850;
const double XMUV3max = 1.200;
const double YMUV3max = 1.200;


const double TDCRES = 0.09765625;


#endif
