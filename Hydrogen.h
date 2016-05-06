#ifndef __HYDROGEN_H__
#define __HYDROGEN_H__


#include <iostream>


namespace Hydrogen{
  const long double m_omega_eg_k = 4161.166504;   // energy separation [kayser]
  // Energy digram in kayser



  const int nOmegaB = 37;
  const int nOmegaC = 14;

  
  // Energy level diagrams
  const long double m_omegaB[nOmegaB] = { 
    90203.5234375000E0 ,
    91521.8906250000E0 ,
    92803.4062500000E0 ,
    94050.0625000000E0 ,
    95263.0781250000E0 ,
    96443.1250000000E0 ,
    97590.6171875000E0 ,
    98705.9687500000E0 ,
    99789.8281250000E0 ,
    100843.359375000E0 ,
    101868.484375000E0 ,
    102868.179687500E0 ,
    103846.679687500E0 ,
    104809.765625000E0 ,
    105765.046875000E0 ,
    106722.171875000E0 ,
    107693.132812500E0 ,
    108692.484375000E0 ,
    109737.640625000E0 ,
    110849.093750000E0 ,
    112050.734375000E0 ,
    113370.023437500E0 ,
    114838.343750000E0 ,
    116491.187500000E0 ,
    118368.453125000E0 ,
    120514.703125000E0 ,
    122979.406250000E0 ,
    125817.203125000E0 ,
    129088.164062500E0 ,
    132858.078125000E0 ,
    137198.656250000E0 ,
    142187.828125000E0 ,
    147910.015625000E0 ,
    154456.328125000E0 ,
    161924.906250000E0 ,
    170421.125000000E0 ,
    180057.843750000E0
  };



const long double m_omegaC[nOmegaC] ={
  99120.1718750000E0 ,
  101427.062500000E0 ,
  103600.281250000E0 ,
  105642.226562500E0 ,
  107554.296875000E0 ,
  109336.890625000E0 ,
  110989.421875000E0 ,
  112510.296875000E0 ,
  113896.929687500E0 ,
  115145.742187500E0 ,
  116252.164062500E0 ,
  117210.609375000E0 ,
  118014.515625000E0 ,
  118656.328125000E0
};

const long double m_uaB[nOmegaB] ={
  1.477046204141564E-061 ,
  4.990470742951365E-061 ,
  9.826115176762043E-061 ,
  1.471993748118731E-060 ,
  1.863141506428632E-060 ,
  2.102900554999763E-060 ,
  2.185679886948450E-060 ,
  2.136207702466737E-060 ,
  1.994430509305663E-060 ,
  1.797590038311674E-060 ,
  1.575841385888343E-060 ,
  1.352712362680451E-060 ,
  1.142465657641083E-060 ,
  9.528397994182980E-061 ,
  7.868613734399500E-061 ,
  6.452789856636180E-061 ,
  5.262919106634739E-061 ,
  4.275406767916083E-061 ,
  3.464809430818794E-061 ,
  2.803120497125353E-061 ,
  2.266178957519624E-061 ,
  1.831367832606562E-061 ,
  1.479608139766074E-061 ,
  1.195874374871601E-061 ,
  9.663177637794816E-062 ,
  7.808872673513472E-062 ,
  6.295726385470742E-062 ,
  5.051510616526656E-062 ,
  4.034993191706247E-062 ,
  3.225222854411815E-062 ,
  2.591922326712091E-062 ,
  2.053817864028401E-062 ,
  1.579169536839474E-062 ,
  1.160864853646887E-062 ,
  7.872553765167933E-063 ,
  4.004801694609017E-063 ,
  5.449997066086150E-064
};
  
const long double m_ubB[nOmegaB] ={
  1.177933816129068E-060 ,
  2.797023543226959E-060 ,
  3.763842173078928E-060 ,
  3.707261350406200E-060 ,
  2.916800966845596E-060 ,
  1.881209714639465E-060 ,
  9.725306997863589E-061 ,
  3.631501885349147E-061 ,
  6.393012017133262E-062 ,
  2.227981789837215E-063 ,
  8.500378196206511E-062 ,
  2.318126883598048E-061 ,
  3.865186905940859E-061 ,
  5.171408588004794E-061 ,
  6.098599360982753E-061 ,
  6.626942532176997E-061 ,
  6.803256305192789E-061 ,
  6.700494207546838E-061 ,
  6.396594747505346E-061 ,
  5.961048304900244E-061 ,
  5.451631198584952E-061 ,
  4.910542973615381E-061 ,
  4.368630781838950E-061 ,
  3.845833729588789E-061 ,
  3.356501051378960E-061 ,
  2.905478822850629E-061 ,
  2.492716284971985E-061 ,
  2.114936314659782E-061 ,
  1.775444442674687E-061 ,
  1.484414726526084E-061 ,
  1.241518337524164E-061 ,
  1.019491764342495E-061 ,
  8.088097893987776E-062 ,
  6.098400311681315E-062 ,
  4.219797623417943E-062 ,
  2.174777961068070E-062 ,
  2.968656691661420E-064
};

const long double m_saB[nOmegaB] ={
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0
};

const long double m_sbB[nOmegaB] = {
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0
};

const long double m_uaC[nOmegaC] = {
  3.788199473966809E-060 ,
  5.665053589688569E-060 ,
  5.316267718941628E-060 ,
  4.085982251657287E-060 ,
  2.841317616876678E-060 ,
  1.874391735996285E-060 ,
  1.208245719704907E-060 ,
  7.705367466950005E-061 ,
  4.916000418386709E-061 ,
  3.145878382677154E-061 ,
  2.019388198584167E-061 ,
  1.284774587260301E-061 ,
  7.827253509524888E-062 ,
  3.716283492687877E-062
};

const long double m_ubC[nOmegaC] = {
  1.089765794878541E-059 ,
  4.536266994221692E-060 ,
  3.915216285046766E-061 ,
  2.067452335366197E-061 ,
  1.116072515853721E-060 ,
  1.714357150633316E-060 ,
  1.822218031056469E-060 ,
  1.630664298473603E-060 ,
  1.326862938045571E-060 ,
  1.018419321616016E-060 ,
  7.501186773071980E-061 ,
  5.291438869534894E-061 ,
  3.472237700586684E-061 ,
  1.727765664496428E-061
};

const long double m_saC[nOmegaC] = {
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0 ,
  1.00000000000000E0 ,
  -1.00000000000000E0
};

 const long double m_sbC[nOmegaC]= {
   1.00000000000000E0 ,
   -1.00000000000000E0 ,
   1.00000000000000E0 ,
   1.00000000000000E0 ,
   -1.00000000000000E0 ,
   1.00000000000000E0 ,
   -1.00000000000000E0 ,
   1.00000000000000E0 ,
   -1.00000000000000E0 ,
   1.00000000000000E0 ,
   -1.00000000000000E0 ,
   1.00000000000000E0 ,
   -1.00000000000000E0 ,
   1.00000000000000E0
 };  
};


#endif


