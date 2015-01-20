<h1>SoundingFilter</h1>
A function to filter high-resolution sounding data in a manner similar to routinely-available sounding data.  It's not perfect, but it's close.

*Authors*:

Greg Blumberg

Tim Supinie

<h1>Procedures</h1>
<h2>Mandatory Levels</h2>

<h3>Mandatory pressure levels</h3>
<ul>
<li>1000 mb
<li>925 mb
<li>850 mb
<li>700 mb
<li>500 mb
<li>400 mb
<li>300 mb
<li>250 mb
<li>200 mb
<li>150 mb
<li>100 mb
<li>70 mb 
<li>50 mb
<li>30 mb
<li>20 mb
<li>10 mb
</ul>

<h3>Tropopause level</h3>
When more than one tropopause is observed, each shall be reported.

The tropopause is defined as the lowest level at which the lapse rate decreases to 2˚C km-1 or less, provided that the average lapse rate between this level and all higher levels within 2 km does not exceed 2˚C km-1. Occasionally, a second tropopause may be found if the lapse rate above the first tropopause exceeds 3˚C km-1. 

<h3>Level of maximum wind</h3>

For coding purposes, a maximum wind level:

(a) Shall be determined by consideration of the list of significant levels for wind speed, as
obtained by means of the relevant recommended or equivalent national method and not by consideration of the original wind-speed curve;

(b) Shall be located above the 500-hPa isobaric surface and shall correspond to a speed
of more than 30 metres per second.
*Note:* A maximum wind level is defined as a level at which the wind speed is greater than that observed immediately above and below that level.

Whenever more than one maximum wind level exists, these levels shall be reported as follows:

(a) The level of greatest maximum wind speed shall be transmitted first;

(b) The other levels shall be classified in descending order of speed and be transmitted
only if their speed exceeds those of the two adjacent minimals by at least 10 metres
per second;

(c) The levels of maximum wind with the same speed shall be encoded successively,
beginning with the lowest one;

(d) Furthermore, the highest level attained by the sounding shall be transmitted, provided:

(i) It satisfies the criteria set forth in Regulation 32.2.3.1 above [the one listed above this section];

(ii) It constitutes the level of the greatest speed of the whole sounding.

<h2>Significant Levels</h2>

The reported significant data alone shall make it possible to reconstruct the wind profile with sufficient accuracy for practical use. Care shall be taken that:

(a) The direction and speed curves (in function of the log of pressure or altitude) can be
reproduced with their prominent characteristics;

(b) These curves can be reproduced with an accuracy of at least 10° for direction and five
metres per second for speed;

(c) The number of significant levels is kept strictly to a necessary minimum.

*Note:* To satisfy these criteria, the following method of successive approximations is recommended, but other methods of attaining equivalent results may suit some national practices better and may be used:

(1) The surface level and the highest level attained by the sounding constitute the first and the last significant levels.
The deviation from the linearly interpolated values between these two levels is then considered. If no direction deviates by more than 10° and no speed by more than five metres per second, no other significant level need be reported. Whenever one parameter deviates by more than the limit specified in paragraph (b) above, the level of greatest deviation becomes a supplementary significant level for both parameters.

(2) The additional significant levels so introduced divide the sounding into two layers. In each separate layer, the deviations from the linearly interpolated values between the base and the top are then considered. The process used in paragraph (1) above is repeated and yields other significant levels. These additional levels in turn modify the layer distribution, and the method is applied again until any level is approximated to the above-mentioned specified values.

For the purpose of computational work, it should be noted that the values derived from a PILOT report present two different resolutions:

(a) Winds at significant levels are reported to the resolution of 5° in direction and one metre per
second in speed;

(b) Any interpolated wind at a level between two significant levels is implicitly reported to the
resolution of ± 10° in direction and ± 5 metres per second in speed.

<h2>Mandatory Significant Levels</h2>

The following shall be included as “mandatory significant levels”:

(a) Surface level and the highest level of the sounding, or aircraft reference level and ter-
mination level for descent soundings;

(b) A level between 110 and 100 hPa;

(c) Bases and tops of inversions and isothermal layers which are at least 20 hPa thick,
provided that the base of the layer occurs below the 300-hPa level or the first
tropopause, whichever is the higher;

(d) Bases and tops of inversion layers which are characterized by a change in temperature of at least 2.5°C or a change in relative humidity of at least 20 per cent, provided that the base of the layer occurs below the 300-hPa level or the first tropopause, whichever is the higher.

*Note:* The inversion layers of (c) and (d) may be comprised of several thinner inversion layers separated by thin layers of temperature lapse. To allow for this situation, the tops of the inversion layers of (c) and (d) shall each be at a level such that no further inversion layers, whether thick or thin, shall occur for at least 20 hPa above the level.

*Note:* Criteria for determining significant levels with respect to changes in wind speed and direction are given in Regulation 32.3.1.

<h2>Additional Levels</h2>

The following shall be included as “additional levels”. They shall be selected in the order given, thereby giving priority to representing the temperature profile. As far as possible, these additional levels shall be the actual levels at which prominent changes in the lapse rate of air temperature occur:

(a) Levels which are necessary to ensure that the temperature obtained by linear interpolation (on a T-log P or essentially similar diagram) between adjacent significant levels shall not depart from the observed temperature by more than 1°C below the first significant level reported above the 300-hPa level or the first tropopause, whichever level is the lower, or by more than 2°C thereafter;

(b) Levels which are necessary to ensure that the relative humidity obtained by linear interpolation between adjacent significant levels shall not depart by more than 15 per cent from the observed values. (The criterion of 15 per cent refers to an amount of relative humidity and NOT to the percentage of the observed value, e.g. if an observed value is 50 per cent, the interpolated value shall lie between 35 per cent and 65 per cent.);

(c) Levels which are necessary to limit the interpolation error on diagrams other than T-log P. These levels shall be such that the pressure at one significant level divided by the pressure of the preceding significant level shall exceed 0.6 for levels up to the first tropopause and shall be determined by use of the method for selecting additional levels but with application of tighter criteria.

[With respect to extrapolating sounding data]

Whenever it is desired to extrapolate a sounding for the computation of the geopotential at
a standard isobaric surface, the following rules shall apply:

(a) Extrapolation is permissible if, and only if, the pressure difference between the minimum pressure of the sounding and the isobaric surface for which the extrapolated value is being computed does not exceed one quarter of the pressure at which the extrapolated value is desired, provided the extrapolation does not extend through a pressure interval exceeding 25 hPa;
(b) For the purpose of geopotential calculation, and for this purpose only, the sounding will be extrapolated, using two points only of the sounding curve on a T-log p diagram, namely that at the minimum pressure reached by the sounding and that at the press- ure given by the sum of this minimum pressure and the pressure difference, mentioned in (a) above.

Paraphrased and parsed from: http://www.wmo.int/pages/prog/www/WMOCodes/Manual/WMO306_Vol-I-1-PartA.pdf
