package io.github.rocsg.fijirelax.mrialgo;

/** This enum define cases of noise management. NOTHING means no noise is estimated, OFFSET means noise is estimated as an additional BIAS, and RICE means full Rice Noise estimation 
 * See the paper Fernandez et al. 2022 (in prep) for more information
 *
 */
public enum NoiseManagement {
	NOTHING,
	OFFSET,
	RICE
}