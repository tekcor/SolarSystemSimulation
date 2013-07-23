//transforms
float4x4 tW: WORLD;        //the models world matrix
float4x4 tV: VIEW;         //view matrix as set via Renderer (EX9)
float4x4 tP: PROJECTION;   //projection matrix as set via Renderer (EX9)
float4x4 tWVP: WORLDVIEWPROJECTION;
float4x4 tWIT : WORLDINVERSETRANSPOSE;
float4x4 tVI : VIEWINVERSE;
float4x4 tWI: WORLDINVERSE;
//material properties

texture tex0;
sampler s0=sampler_state{Texture=(tex0);MipFilter=LINEAR;MinFilter=LINEAR;MagFilter=LINEAR;AddressU=CLAMP;AddressV=CLAMP;};

float4x4 tTex: TEXTUREMATRIX <string uiname="Texture Transform";>;

struct vs2ps
{
    float4 Pos : POSITION;
	float3 Norm : NORMAL0;
	float4 TexCd : TEXCOORD0;
	float3 View : TEXCOORD1;
	float3 PosW : TEXCOORD2;
	float4 c0 : COLOR0;
	float4 c1 : COLOR1;
};
float4 Color:COLOR =1;
//uniform float3 v3CameraPos;		// The camera's current position
float3 lPos={1,0,0};
//#define v3LightPos (normalize(lPos))
//uniform float3 v3LightPos;		// The direction vector to the light source
uniform float3 v3InvWavelength;	// 1 / pow(wavelength, 4) for the red, green, and blue channels
//uniform float fCameraHeight;	// The camera's current height
//#define fCameraHeight2 (fCameraHeight*fCameraHeight)
//uniform float fCameraHeight2;	// fCameraHeight^2
uniform float fOuterRadius;		// The outer (atmosphere) radius
#define fOuterRadius2 (fOuterRadius*fOuterRadius)
//uniform float fOuterRadius2;	// fOuterRadius^2
uniform float fInnerRadius;		// The inner (planetary) radius
#define fInnerRadius2 (fInnerRadius*fInnerRadius)
//uniform float fInnerRadius2;	// fInnerRadius^2
uniform float fESun=15;
uniform float fKr=.0025;
uniform float fKm=.0015;
#define fKrESun (fKr*fESun)
#define fKmESun (fKm*fESun)
#define fKr4PI (fKr*4*acos(-1))
#define fKm4PI (fKm*4*acos(-1))
//uniform float fKrESun;			// Kr * ESun
//uniform float fKmESun;			// Km * ESun
//uniform float fKr4PI;			// Kr * 4 * PI
//uniform float fKm4PI;			// Km * 4 * PI
#define fScale (1.0/(fOuterRadius - fInnerRadius))
//uniform float fScale;			// 1 / (fOuterRadius - fInnerRadius)
uniform float fScaleDepth;		// The scale depth (i.e. the altitude at which the atmosphere's average density is found)
#define fScaleOverScaleDepth (fScale/fScaleDepth)
//uniform float fScaleOverScaleDepth;	// fScale / fScaleDepth

const int nSamples = 2;
//const float fSamples = 2.0;

uniform float g;
#define g2 (g*g)
//uniform float g2;
float GammaFactor <float uimin=0.0; float uimax=1.0;> = 1;
float SunGlow <float uimin=0.0;> =1;
float scale(float fCos)
{
	float x = 1.0 - fCos;
	return fScaleDepth * exp(-0.00287 + x*(0.459 + x*(3.83 + x*(-6.80 + x*5.25))));
}
//*
vs2ps VS_Ground(
    float4 Pos : POSITION,
    float4 TexCd : TEXCOORD0,
	float3 Norm : NORMAL0)
{
	vs2ps Out = (vs2ps)0;
	Out.Pos = mul(Pos, tWVP);
	Out.TexCd=TexCd;
	float4 PosW=mul(Pos,tW);
	
	float3 v3CameraPos=mul(float4(0,0,0,1),tVI).xyz;
	float3 v3Pos = PosW.xyz;
	float3 v3Ray = v3Pos - v3CameraPos;
	
	//Recalc vectors and positions for world transform
	float3 v3Direction;
	v3Pos=mul(float4(v3Pos.xyz,1),tWI);
	v3CameraPos=mul(float4(v3CameraPos.xyz,1),tWI);
	v3Ray = v3Pos - v3CameraPos;
	v3Direction=v3Pos - v3CameraPos;
	float fCameraHeight=distance(v3CameraPos,0);
	float fCameraHeight2=fCameraHeight*fCameraHeight;
	float3 v3LightPos=normalize(mul(float4(lPos,1),tWI).xyz);
	
	float fFar = length(v3Ray);
	v3Ray /= fFar;
	
	// Calculate the closest intersection of the ray with the outer atmosphere (which is the near point of the ray passing through the atmosphere)
	float B = 2.0 * dot(v3CameraPos, v3Ray);
	float C = fCameraHeight2 - fOuterRadius2;
	float fDet = max(0.0, B*B - 4.0 * C);
	float fNear = 0.5 * (-B - sqrt(fDet));

	// Calculate the ray's starting position, then calculate its scattering offset
	float3 v3Start = v3CameraPos + v3Ray * fNear;
	fFar -= fNear;
	
	if(fCameraHeight<fOuterRadius){
		v3Start = v3CameraPos;
		fFar += fNear;
	}
	// Calculate the ray's starting position, then calculate its scattering offset
	//float3 v3Start = v3CameraPos;
	float fDepth = exp((fInnerRadius - fCameraHeight) / fScaleDepth);
	float fCameraAngle = dot(-v3Ray, v3Pos) / length(v3Pos);
	float fLightAngle = dot(v3LightPos, v3Pos) / length(v3Pos);
	float fCameraScale = scale(fCameraAngle);
	float fLightScale = scale(fLightAngle);
	float fCameraOffset = fDepth*fCameraScale;
	float fTemp = (fLightScale + fCameraScale);

	// Initialize the scattering loop variables
	float fSampleLength = fFar / float(nSamples);
	float fScaledLength = fSampleLength * fScale;
	float3 v3SampleRay = v3Ray * fSampleLength;
	float3 v3SamplePoint = v3Start + v3SampleRay * 0.5;

	// Now loop through the sample rays
	float3 v3FrontColor = float3(0.0, 0.0, 0.0);
	float3 v3Attenuate;
	for(int i=0; i<nSamples; i++)
	{
		float fHeight = length(v3SamplePoint);
		float fDepth = exp(fScaleOverScaleDepth * (fInnerRadius - fHeight));
		float fScatter = fDepth*fTemp - fCameraOffset;
		v3Attenuate = exp(-fScatter * (v3InvWavelength * fKr4PI + fKm4PI));
		v3FrontColor += v3Attenuate * (fDepth * fScaledLength);
		v3SamplePoint += v3SampleRay;
	}

	//gl_FrontColor.rgb = v3FrontColor * (v3InvWavelength * fKrESun + fKmESun);

	// Calculate the attenuation factor for the ground
	//gl_FrontSecondaryColor.rgb = v3Attenuate;

	Out.c0.rgb = v3FrontColor * (v3InvWavelength * fKrESun + fKmESun);
	Out.c1.rgb = v3Attenuate;
	Out.View=PosW.xyz-mul(float4(0,0,0,1),tVI);
	Out.PosW=PosW;
    return Out;
}

vs2ps VS_Sky(
    float4 Pos : POSITION,
    float4 TexCd : TEXCOORD0,
	float3 Norm : NORMAL0)
{
	vs2ps Out = (vs2ps)0;
	Pos.xyz*=fOuterRadius/fInnerRadius;
	Out.Pos = mul(Pos, tWVP);
	Out.TexCd=TexCd;
	float4 PosW=mul(Pos,tW);
	
	// Get the ray from the camera to the vertex and its length (which is the far point of the ray passing through the atmosphere)
	float3 v3CameraPos=mul(float4(0,0,0,1),tVI).xyz;
	float3 v3Pos = PosW.xyz;
	float3 v3Ray = v3Pos - v3CameraPos;
	
	//Recalc vectors and positions for world transform
	float3 v3Direction;
	v3Pos=mul(float4(v3Pos.xyz,1),tWI);
	v3CameraPos=mul(float4(v3CameraPos.xyz,1),tWI);
	v3Ray = v3Pos - v3CameraPos;
	v3Direction=v3Pos - v3CameraPos;
	float fCameraHeight=distance(v3CameraPos,0);
	float fCameraHeight2=fCameraHeight*fCameraHeight;
	float3 v3LightPos=normalize(mul(float4(lPos,1),tWI).xyz);
	
	float fFar = length(v3Ray);
	v3Ray /= fFar;
	
	// Calculate the closest intersection of the ray with the outer atmosphere (which is the near point of the ray passing through the atmosphere)
	float B = 2.0 * dot(v3CameraPos, v3Ray);
	float C = fCameraHeight2 - fOuterRadius2;
	float fDet = max(0.0, B*B - 4.0 * C);
	float fNear = 0.5 * (-B - sqrt(fDet));

	// Calculate the ray's starting position, then calculate its scattering offset
	float3 v3Start = v3CameraPos + v3Ray * fNear;
	fFar -= fNear;
	float fStartAngle = dot(v3Ray, v3Start) / fOuterRadius;
	float fStartDepth = exp(-1.0 / fScaleDepth);
	float fStartOffset = fStartDepth*scale(fStartAngle);
	
	if(fCameraHeight<fOuterRadius){
		v3Start = v3CameraPos;
		float fHeight = length(v3Start);
		float fDepth = exp(fScaleOverScaleDepth * (fInnerRadius - fCameraHeight));
		fStartAngle = dot(v3Ray, v3Start) / fHeight;
		fStartOffset = fDepth*scale(fStartAngle);
		fFar += fNear;
	}

	// Initialize the scattering loop variables
	//gl_FrontColor = vec4(0.0, 0.0, 0.0, 0.0);
	float fSampleLength = fFar / float(nSamples);
	float fScaledLength = fSampleLength * fScale;
	float3 v3SampleRay = v3Ray * fSampleLength;
	float3 v3SamplePoint = v3Start + v3SampleRay * 0.5;

	// Now loop through the sample rays
	float3 v3FrontColor = float3(0.0, 0.0, 0.0);
	for(int i=0; i<nSamples; i++)
	{
		float fHeight = length(v3SamplePoint);
		float fDepth = exp(fScaleOverScaleDepth * (fInnerRadius - fHeight));
		float fLightAngle = dot(v3LightPos, v3SamplePoint) / fHeight;
		float fCameraAngle = dot(v3Ray, v3SamplePoint) / fHeight;
		float fScatter = (fStartOffset + fDepth*(scale(fLightAngle) - scale(fCameraAngle)));
		float3 v3Attenuate = exp(-fScatter * (v3InvWavelength * fKr4PI + fKm4PI));
		v3FrontColor += v3Attenuate * (fDepth * fScaledLength);
		v3SamplePoint += v3SampleRay;
	}

	// Finally, scale the Mie and Rayleigh colors and set up the varying variables for the pixel shader
	//gl_FrontSecondaryColor.rgb = v3FrontColor * fKmESun;
	//gl_FrontColor.rgb = v3FrontColor * (v3InvWavelength * fKrESun);
	Out.c0.rgb = v3FrontColor * (v3InvWavelength * fKrESun);
	Out.c1.rgb = v3FrontColor * fKmESun;
	Out.View=PosW.xyz-mul(float4(0,0,0,1),tVI);
	Out.PosW=PosW;
    return Out;
}

// --------------------------------------------------------------------------------------------------
// PIXELSHADERS:
// --------------------------------------------------------------------------------------------------

float4 PS_Ground(vs2ps In): COLOR
{
	float3 v3Direction=In.View;
	float4 gl_FragColor=1;
	float4 gl_Color=In.c0;
	float4 gl_SecondaryColor=In.c1;
	
	gl_FragColor = gl_Color + 0.25 * gl_SecondaryColor;
	gl_FragColor.a=1;
	
	//calc fCameraHeight
	float3 v3CameraPos=mul(float4(0,0,0,1),tVI).xyz;
	v3CameraPos=mul(float4(v3CameraPos.xyz,1),tWI);
	float fCameraHeight=distance(v3CameraPos,0);
	//Gamma correction inside atmosphere, recalc color
	float Gamma=lerp(1,pow(smoothstep(fInnerRadius,fOuterRadius,fCameraHeight),0.5),GammaFactor);
	gl_FragColor.rgb=normalize(gl_Color)*pow(length(gl_Color.rgb)/sqrt(2),Gamma)*sqrt(2);
	gl_FragColor.rgb+=.25*SunGlow*gl_SecondaryColor;
	
	
	gl_FragColor*=Color;
	//gl_FragColor+=pow(gl_FragColor,.5)*tex2D(s0,In.TexCd.xy);  //texture
    return gl_FragColor;
}

float4 PS_Sky(vs2ps In): COLOR
{
	float3 v3Direction=In.View;
	float4 gl_FragColor=1;
	float4 gl_Color=In.c0;
	float4 gl_SecondaryColor=In.c1;
	float3 v3LightPos=normalize(mul(float4(lPos,1),tWI).xyz);
	
	float fCos = dot(v3LightPos, v3Direction) / length(v3Direction);
	//float fRayleighPhase = 0.75 * (1.0 + fCos*fCos);
	float fRayleighPhase = 0.75 * (1.0 + fCos);
	//float fMiePhase = 1.5 * ((1.0 - g2) / (2.0 + g2)) * (1.0 + fCos*fCos) / pow(1.0 + g2 - 2.0*g*fCos, 1.5);
	float fMiePhase = 1.5 * ((1.0 - g2) / (2.0 + g2)) * (1.0 + fCos*fCos) / pow(1.0 + g2 + 2.0*g*fCos, 1.5);
	
	gl_FragColor = fRayleighPhase * gl_Color + fMiePhase * gl_SecondaryColor;
	gl_FragColor.a=1;
	
	//calc fCameraHeight
	float3 v3CameraPos=mul(float4(0,0,0,1),tVI).xyz;
	v3CameraPos=mul(float4(v3CameraPos.xyz,1),tWI);
	float fCameraHeight=distance(v3CameraPos,0);
	//Gamma correction inside atmosphere, recalc color
	float Gamma=lerp(1,pow(smoothstep(fInnerRadius,fOuterRadius,fCameraHeight),0.5),GammaFactor);
	gl_FragColor.rgb=normalize(gl_Color)*pow(fRayleighPhase,pow(Gamma,.5))*pow(length(gl_Color.rgb)/sqrt(2),Gamma)*sqrt(2);
	gl_FragColor.rgb+=SunGlow*fMiePhase * gl_SecondaryColor;
	
	gl_FragColor*=Color;
	//gl_FragColor.rgb=fRayleighPhase;
    return gl_FragColor;
}

//*/
vs2ps VS_GroundHQ(
    float4 Pos : POSITION,
    float4 TexCd : TEXCOORD0,
	float3 Norm : NORMAL0)
{
	vs2ps Out = (vs2ps)0;
	Out.Pos = mul(Pos, tWVP);
	Out.TexCd=TexCd;
	float4 PosW=mul(Pos,tW);
	
	Out.View=PosW.xyz-mul(float4(0,0,0,1),tVI);
	Out.PosW=PosW;
    return Out;
}
vs2ps VS_SkyHQ(
    float4 Pos : POSITION,
    float4 TexCd : TEXCOORD0,
	float3 Norm : NORMAL0)
{
	vs2ps Out = (vs2ps)0;
	Pos.xyz*=fOuterRadius/fInnerRadius;
	Out.Pos = mul(Pos, tWVP);
	Out.TexCd=TexCd;
	float4 PosW=mul(Pos,tW);
	Out.View=PosW.xyz-mul(float4(0,0,0,1),tVI);
	Out.PosW=PosW;
    return Out;
}

float4 PS_GroundHQ(vs2ps In): COLOR
{
	float3 PosW=In.PosW;
	
	float3 v3CameraPos=mul(float4(0,0,0,1),tVI).xyz;
	float3 v3Pos = PosW.xyz;
	float3 v3Ray = v3Pos - v3CameraPos;
	
	//Recalc vectors and positions for world transform
	float3 v3Direction;
	v3Pos=mul(float4(v3Pos.xyz,1),tWI);
	v3CameraPos=mul(float4(v3CameraPos.xyz,1),tWI);
	v3Ray = v3Pos - v3CameraPos;
	v3Direction=v3Pos - v3CameraPos;
	float fCameraHeight=distance(v3CameraPos,0);
	float fCameraHeight2=fCameraHeight*fCameraHeight;
	float3 v3LightPos=normalize(mul(float4(lPos,1),tWI).xyz);
	
	float fFar = length(v3Ray);
	v3Ray /= fFar;
	
	// Calculate the closest intersection of the ray with the outer atmosphere (which is the near point of the ray passing through the atmosphere)
	float B = 2.0 * dot(v3CameraPos, v3Ray);
	float C = fCameraHeight2 - fOuterRadius2;
	float fDet = max(0.0, B*B - 4.0 * C);
	float fNear = 0.5 * (-B - sqrt(fDet));

	// Calculate the ray's starting position, then calculate its scattering offset
	float3 v3Start = v3CameraPos + v3Ray * fNear;
	fFar -= fNear;
	
	if(fCameraHeight<fOuterRadius){
		v3Start = v3CameraPos;
		fFar += fNear;
	}
	// Calculate the ray's starting position, then calculate its scattering offset
	//float3 v3Start = v3CameraPos;
	float fDepth = exp((fInnerRadius - fCameraHeight) / fScaleDepth);
	float fCameraAngle = dot(-v3Ray, v3Pos) / length(v3Pos);
	float fLightAngle = dot(v3LightPos, v3Pos) / length(v3Pos);
	float fCameraScale = scale(fCameraAngle);
	float fLightScale = scale(fLightAngle);
	float fCameraOffset = fDepth*fCameraScale;
	float fTemp = (fLightScale + fCameraScale);

	// Initialize the scattering loop variables
	float fSampleLength = fFar / float(nSamples);
	float fScaledLength = fSampleLength * fScale;
	float3 v3SampleRay = v3Ray * fSampleLength;
	float3 v3SamplePoint = v3Start + v3SampleRay * 0.5;

	// Now loop through the sample rays
	float3 v3FrontColor = float3(0.0, 0.0, 0.0);
	float3 v3Attenuate;
	for(int i=0; i<nSamples; i++)
	{
		float fHeight = length(v3SamplePoint);
		float fDepth = exp(fScaleOverScaleDepth * (fInnerRadius - fHeight));
		float fScatter = fDepth*fTemp - fCameraOffset;
		v3Attenuate = exp(-fScatter * (v3InvWavelength * fKr4PI + fKm4PI));
		v3FrontColor += v3Attenuate * (fDepth * fScaledLength);
		v3SamplePoint += v3SampleRay;
	}

	//gl_FrontColor.rgb = v3FrontColor * (v3InvWavelength * fKrESun + fKmESun);

	// Calculate the attenuation factor for the ground
	//gl_FrontSecondaryColor.rgb = v3Attenuate;

	In.c0.rgb = v3FrontColor * (v3InvWavelength * fKrESun + fKmESun);
	In.c1.rgb = v3Attenuate;
	
	float4 gl_FragColor=1;
	float4 gl_Color=In.c0;
	float4 gl_SecondaryColor=In.c1;
	
	gl_FragColor = gl_Color + 0.25 * gl_SecondaryColor;
	gl_FragColor.a=1;
	
	//Gamma correction inside atmosphere, recalc color
	float Gamma=lerp(1,pow(smoothstep(fInnerRadius,fOuterRadius,fCameraHeight),0.5),GammaFactor);
	gl_FragColor.rgb=normalize(gl_Color)*pow(length(gl_Color.rgb)/sqrt(2),Gamma)*sqrt(2);
	gl_FragColor.rgb+=.25*SunGlow*gl_SecondaryColor;
	
	gl_FragColor*=Color;
	//gl_FragColor+=pow(gl_FragColor,.5)*tex2D(s0,In.TexCd.xy);
    return gl_FragColor;
}
float4 PS_SkyHQ(vs2ps In): COLOR
{
	float3 PosW=In.PosW;
	// Get the ray from the camera to the vertex and its length (which is the far point of the ray passing through the atmosphere)
	float3 v3CameraPos=mul(float4(0,0,0,1),tVI);
	float3 v3Pos = PosW.xyz;
	float3 v3Ray = v3Pos - v3CameraPos;
	
	
	//Recalc vectors and positions for world transform
	float3 v3Direction;
	v3Pos=mul(float4(v3Pos.xyz,1),tWI);
	v3CameraPos=mul(float4(v3CameraPos.xyz,1),tWI);
	v3Ray = v3Pos - v3CameraPos;
	v3Direction=v3Pos - v3CameraPos;
	float fCameraHeight=distance(v3CameraPos,0);
	float fCameraHeight2=fCameraHeight*fCameraHeight;
	float3 v3LightPos=normalize(mul(float4(lPos,1),tWI).xyz);
	
	float fFar = length(v3Ray);
	v3Ray /= fFar;
	
	// Calculate the closest intersection of the ray with the outer atmosphere (which is the near point of the ray passing through the atmosphere)
	float B = 2.0 * dot(v3CameraPos, v3Ray);
	float C = fCameraHeight2 - fOuterRadius2;
	float fDet = max(0.0, B*B - 4.0 * C);
	float fNear = 0.5 * (-B - sqrt(fDet));

	// Calculate the ray's starting position, then calculate its scattering offset
	float3 v3Start = v3CameraPos + v3Ray * fNear;
	fFar -= fNear;
	float fStartAngle = dot(v3Ray, v3Start) / fOuterRadius;
	float fStartDepth = exp(-1.0 / fScaleDepth);
	float fStartOffset = fStartDepth*scale(fStartAngle);
	
	if(fCameraHeight<fOuterRadius){
		v3Start = v3CameraPos;
		float fHeight = length(v3Start);
		float fDepth = exp(fScaleOverScaleDepth * (fInnerRadius - fCameraHeight));
		fStartAngle = dot(v3Ray, v3Start) / fHeight;
		fStartOffset = fDepth*scale(fStartAngle);
		fFar += fNear;
	}

	// Initialize the scattering loop variables
	//gl_FrontColor = vec4(0.0, 0.0, 0.0, 0.0);
	float fSampleLength = fFar / (float)nSamples;
	float fScaledLength = fSampleLength * fScale;
	float3 v3SampleRay = v3Ray * fSampleLength;
	float3 v3SamplePoint = v3Start + v3SampleRay * 0.5;

	// Now loop through the sample rays
	float3 v3FrontColor = float3(0.0, 0.0, 0.0);
	for(int i=0; i<nSamples; i++)
	{
		float fHeight = length(v3SamplePoint);
		float fDepth = exp(fScaleOverScaleDepth * (fInnerRadius - fHeight));
		float fLightAngle = dot(v3LightPos, v3SamplePoint) / fHeight;
		float fCameraAngle = dot(v3Ray, v3SamplePoint) / fHeight;
		float fScatter = (fStartOffset + fDepth*(scale(fLightAngle) - scale(fCameraAngle)));
		float3 v3Attenuate = exp(-fScatter * (v3InvWavelength * fKr4PI + fKm4PI));
		v3FrontColor += v3Attenuate * (fDepth * fScaledLength);
		v3SamplePoint += v3SampleRay;
	}

	// Finally, scale the Mie and Rayleigh colors and set up the varying variables for the pixel shader
	//gl_FrontSecondaryColor.rgb = v3FrontColor * fKmESun;
	//gl_FrontColor.rgb = v3FrontColor * (v3InvWavelength * fKrESun);
	In.c0.rgb = v3FrontColor * (v3InvWavelength * fKrESun);
	In.c1.rgb = v3FrontColor * fKmESun;
	
	float4 gl_FragColor=1;
	float4 gl_Color=In.c0;
	float4 gl_SecondaryColor=In.c1;
	float fCos = dot(v3LightPos, v3Direction) / length(v3Direction);
	//float fRayleighPhase = 0.75 * (1.0 + fCos*fCos);
	float fRayleighPhase = 0.75 * (1.0 + fCos);
	//float fMiePhase = 1.5 * ((1.0 - g2) / (2.0 + g2)) * (1.0 + fCos*fCos) / pow(1.0 + g2 - 2.0*g*fCos, 1.5);
	
	float fMiePhase = 1.5 * ((1.0 - g2) / (2.0 + g2)) * (1.0 + fCos*fCos) / pow(1.0 + g2 + 2.0*g*fCos, 1.5);
	gl_FragColor = fRayleighPhase * gl_Color + fMiePhase * gl_SecondaryColor;
	gl_FragColor.a=1;
	
	
	//Gamma correction inside atmosphere, recalc color
	float Gamma=lerp(1,pow(smoothstep(fInnerRadius,fOuterRadius,fCameraHeight),0.5),GammaFactor);
	gl_FragColor.rgb=normalize(gl_Color)*pow(fRayleighPhase,pow(Gamma,.5))*pow(length(gl_Color.rgb)/sqrt(2),Gamma)*sqrt(2);
	gl_FragColor.rgb+=SunGlow*fMiePhase * gl_SecondaryColor;
	
	gl_FragColor*=Color;
	
    return gl_FragColor;
}

// --------------------------------------------------------------------------------------------------
// TECHNIQUES:
// --------------------------------------------------------------------------------------------------

technique _Vertex
{
	

	pass P1
    {
        CullMode=CCW;
    	AlphaBlendEnable=FALSE;
    	ZWriteEnable=TRUE;
        VertexShader = compile vs_3_0 VS_Ground();
        PixelShader = compile ps_3_0 PS_Ground();
    }
    pass P3
    {
    	CullMode=CW;
    	AlphaBlendEnable=TRUE;
		SrcBlend=ONE;
		DestBlend=ONE;
    	ZWriteEnable=FALSE;
        VertexShader = compile vs_3_0 VS_Sky();
        PixelShader = compile ps_3_0 PS_Sky();
    }
}

technique _Pixel
{
	

	pass P1
    {
        CullMode=CCW;
    	AlphaBlendEnable=FALSE;
    	ZWriteEnable=TRUE;
        VertexShader = compile vs_3_0 VS_GroundHQ();
        PixelShader = compile ps_3_0 PS_GroundHQ();
    }
    pass P3
    {
    	CullMode=CW;
    	AlphaBlendEnable=TRUE;
		SrcBlend=ONE;
		DestBlend=ONE;
    	ZWriteEnable=FALSE;
        VertexShader = compile vs_3_0 VS_SkyHQ();
        PixelShader = compile ps_3_0 PS_SkyHQ();
    }
}