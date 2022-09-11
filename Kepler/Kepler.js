//loome eraldi massiivid koefitsentide jaoks
var a=new Array(6);
var b=new Array(6);
var ch=new Array(6);
var ct=new Array(6);
var k=new Array(6);
var y=new Array(4);
var ynew=new Array(4);
var te=new Array(4);
var yy=new Array(4);
var yp=new Array(4);

var graafikute_kõrgus=300;
var plot_offset_X=100;
var plot_offset_Y=graafikute_kõrgus/2;

var n=4; // võrrandite arv diff võrrandi süsteemis
var e=0.025; // suvaline error
var h=100;
var t=0.0;
var nt=600000.0;
var it=0;

// Keskkond
var g_acceleration = 9.8; // m/s^2
var G_const = 6.67430e-11;// kg^-2 N*m^2


var model_num=["Maa","Mars","Veenus", "Jupiter"];
var planet_radius=[6371.0, 3389.5, 6051.8, 69911.0]; // in kilometers
var planet_mass=[5.97219e24, 6.39e23, 4.867e24, 1.898e27]; //in kilogramm
//var planet_mass=[5.97219*pow(10,24), 6.39*pow(10,23), 4.867*pow(10,24), 1.898*pow(10,27)]; //in kilogramm
//var planet_mass=[5.97219*10**24, 6.39*10**23, 4.867*10**24, 1.898*10**27];
var gravity_coef = 1.0;
var mass_default = 5.97219e24;
var radius_default = 6371.0;

// Keha
var m_body=100.0; // kilogramm
var r_body=0.001; // kilomeeter
var c_body= 0.47; // 

// k_{i,n} masiiv, i=k_{i} ja n= võrrandi number
for (var i=0; i<6; i++){
k[i]= new Array(n);}

for (var i=0; i<6;i++){
b[i]=new Array(5);}


function coefficients(){

//Fehlebrgi A koefitsendid
a[0]=0;
a[1]=1.0/4.0;
a[2]=3.0/8.0 ;
a[3]=12.0/13.0;
a[4]=1;
a[5]=1.0/2.0;

 // Alguses oli massiv b initsialiseeritud kui 1 veeruline maatriks, nüüd on tal 6 rida ja 5 veergu
// Täidame B massiivi ära Fehlbergi koefitsentidega
b[0][0]=0;
b[0][1]=1.0/4.0;
b[0][2]=3.0/32.0;
b[0][3]=1932.0/2197.0;
b[0][4]=439.0/216.0;
b[0][5]=-8.0/27.0;

b[1][0]=0;
b[1][1]=0;
b[1][2]=9.0/32.0;
b[1][3]=-7200.0/2197.0;
b[1][4]=-8.0;
b[1][5]=2.0;

b[2][0]=0;
b[2][1]=0;
b[2][2]=0;
b[2][3]=7296.0/2197.0;
b[2][4]=3680.0/513.0;
b[2][5]=-3544.0/2565.0;

b[3][0]=0;
b[3][1]=0;
b[3][2]=0;
b[3][3]=0;
b[3][4]=-845.0/4104.0;
b[3][5]=1859.0/4104.0;

b[4][0]=0;
b[4][1]=0;
b[4][2]=0;
b[4][3]=0;
b[4][4]=0;
b[4][5]=-11.0/40.0;

//Fehlbergi CH koefitsendid
ch[0]=16.0/135.0;
ch[1]=0.0;
ch[2]=6656.0/12825.0;
ch[3]=28561.0/56430.0;
ch[4]=-9.0/50.0;
ch[5]=2.0/55.0;
//Fehlbergi CT koefitsendid
ct[0]=1.0/360.0;
ct[1]=0.0;
ct[2]=-128.0/4275.0;
ct[3]=-2197.0/75240;
ct[4]=1.0/50.0;
ct[5]=2.0/55.0;

}


// Differential equations
function diff(tf,yf,ypf) // Neid argumente ju funktsiooni sees ei ole?
{
yp[0]=yy[2];
yp[1]=yy[3];
nimetaja = Math.sqrt(yy[0]*yy[0] + yy[1]*yy[1])*Math.sqrt(yy[0]*yy[0] + yy[1]*yy[1])*Math.sqrt(yy[0]*yy[0] + yy[1]*yy[1]);
yp[2]=-(G_const*mass*1e-9*yy[0])/nimetaja;
yp[3]=-(G_const*(mass*1e-9)*yy[1])/nimetaja;
}


// ------------------------------------------- beginning of integration algorithm -----------------------------
function integration()
{
var neq,i,ynew=new Array(4),te=new Array(4);

 //document.write(neq,"<br>")
  
 //k0
 for (i=0;i<n;i++)
   {yy[i]=y[i];}
 diff();
 for (neq=0;neq<n;neq++)
   {k[0][neq]=h*yp[neq];}
//console.log(yy[0],yy[1])
  
 //k1
 for (i=0;i<n;i++)
   {yy[i]=y[i]+b[0][1]*k[0][i];}
 diff();
 for (neq=0;neq<n;neq++)
   {k[1][neq]=h*yp[neq];}
  
 //k2
 for (i=0;i<n;i++)
   {yy[i]=y[i]+b[0][2]*k[0][i]+b[1][2]*k[1][i];}
 diff();
 for (neq=0;neq<n;neq++)
   {k[2][neq]=h*yp[neq];}
  
 //k3
 for (i=0;i<n;i++)
 {
 yy[i]=y[i]+b[0][3]*k[0][i]+b[1][3]*k[1][i]+b[2][3]*k[2][i];
 }
 diff();
 for (neq=0;neq<n;neq++)
 {k[3][neq]=h*yp[neq];}
  
 //k4
 for (i=0;i<n;i++)
 {
 yy[i]=y[i]+b[0][4]*k[0][i]+b[1][4]*k[1][i]+b[2][4]*k[2][i]+b[3][4]*k[3][i];
 }
 diff();
 for (neq=0;neq<n;neq++)
 {k[4][neq]=h*yp[neq];}
  
 //k5
 for (i=0;i<n;i++)
 {
 yy[i]=y[i]+b[0][5]*k[0][i]+b[1][5]*k[1][i]+b[2][5]*k[2][i]+b[3][5]*k[3][i]+b[4][5]*k[4][i];
 }
 diff();
 for (neq=0;neq<n;neq++)
 {k[5][neq]=h*yp[neq];}
  
 // lopp liitmine
 for (neq=0;neq<n;neq++)
 {ynew[neq]=y[neq]+ch[0]*k[0][neq]+ch[1]*k[1][neq]+ch[2]*k[2][neq]+ch[3]*k[3][neq]+ch[4]*k[4][neq]+ch[5]*k[5][neq];}
 te[neq]=Math.abs(ct[0]*k[neq][0]+ct[1]*k[neq][1]+ct[2]*k[neq][2]+ct[3]*k[neq][3]+ct[4]*k[neq][4]+ct[5]*k[neq][5]);

for (neq=0; neq<n;neq++){
h_uus=0.9*h*pow((e/te[neq]), 1.0/4.0); 
if (te[neq]>e ) {
h=h_uus;
neq=neq-1; // Kordame iteratsiooni järgmisel iteratsioonil
} else if (te[neq]<= e){
h=h_uus; // iteratsiooni ei korda, lähme edasi.
}
}
  
 t=t+h;
 y[0]=ynew[0];
 y[1]=ynew[1];
 y[2]=ynew[2];
 y[3]=ynew[3];
 //console.log(t,y[0],y[1]);
}

// ------------------------------------------- end of integration algorithm -----------------------------

var toggle=true;
function toggle_button() {
  if (toggle==true){
    loop();
    toggle= !toggle;
    STOP_CONTINUE_button.html("Pause");
  }
  else if (toggle==false){ 
        noLoop();
        toggle= !toggle;
        STOP_CONTINUE_button.html("Continue");
   }
}


// --------------------------------------------- INITIALIZING THE MODEL ---------------------------------------------------------
var mudeli_init_loendur=0;
function init_model(num){
if (mudeli_init_loendur>0) {
  plot_scale=1.0;
  anim_scale=1.0;

  aja_massiiv=[];
  Vx_massiiv=[];
  Vy_massiiv=[];
  Ax_massiiv=[];
  Ay_massiiv=[];
  Epot_massiiv=[];
  Ekin_massiiv=[];
  Etot_massiiv=[];  
  h_massiiv=[]; //planeedi tsentrist moodulini  
  V_mod_massiiv=[];
  A_mod_massiiv=[];

  model_name=model_num[num];
  radius=planet_radius[num];
  mass=planet_mass[num];

  t=0.0;
  nt=60000.0;
  it=0.0;

  y[0]=radius+100;
  y[1]=0;
  y[2]=0;
  y[3]=float(VELOCITY_INPUT.value());
  background(0);
} else {
    plot_scale=1.0;
    anim_scale=1.0;

    aja_massiiv=[];
    Vx_massiiv=[];
    Vy_massiiv=[];
    Ax_massiiv=[];
    Ay_massiiv=[];
    Epot_massiiv=[];
    Ekin_massiiv=[];
    Etot_massiiv=[];  
    h_massiiv=[]; //planeedi tsentrist moodulini  
    V_mod_massiiv=[];
    A_mod_massiiv=[];

    model_name=model_num[num];
    radius=planet_radius[num];
    mass=planet_mass[num];

    t=0.0;
    nt=60000.0;
    it=0.0;

    y[0]=radius+100;
    y[1]=0;
    y[2]=0;
    y[3]=Math.sqrt((G_const*mass*1e-9)/(radius+100)); // See enamvähem töötab ja annab numbriliselt õige väärtuse
    //y[3]=9;
    background(0);
    VELOCITY_INPUT=createInput();
    VELOCITY_INPUT.position(width+50,height-110);
    VELOCITY_INPUT.size(70,30);
    VELOCITY_INPUT_TEXT=createP("Sisesta esialgne kiirus: ")
    VELOCITY_INPUT_TEXT.position(width+20,height-150);
}
mudeli_init_loendur=mudeli_init_loendur+1;
}


function setup() {
  
aja_massiiv=[];
Vx_massiiv=[];
Vy_massiiv=[];
Ax_massiiv=[];
Ay_massiiv=[];
Epot_massiiv=[];
Ekin_massiiv=[];
Etot_massiiv=[];  
h_massiiv=[]; //planeedi tsentrist moodulini  
V_mod_massiiv=[];
A_mod_massiiv=[];
  
createCanvas(800,500+graafikute_kõrgus);
background(0);
coefficients();
init_model(0);
model_name=model_num[0];

  
  
STOP_CONTINUE_button = createButton("Start simulation");
STOP_CONTINUE_button.size(100,50);
STOP_CONTINUE_button.position(width+50, (height-graafikute_kõrgus)*0+30);
STOP_CONTINUE_button.mousePressed(toggle_button);

EARTH_button = createButton("Maa");
EARTH_button.size(50,50);
EARTH_button.position(width+50, (height-graafikute_kõrgus)*0+100);
EARTH_button.mousePressed(function(){init_model(0);});

MARS_button = createButton("Mars");
MARS_button.size(100,50);
MARS_button.position(width+50, (height-graafikute_kõrgus)*0+160);
MARS_button.mousePressed(function(){init_model(1);});

VENUS_button = createButton("Veenus");
VENUS_button.size(100,50);
VENUS_button.position(width+50, (height-graafikute_kõrgus)*0+220);
VENUS_button.mousePressed(function(){init_model(2);});

// MASS_SLIDER = createSlider(0.1e24, 10e24, 5.97219e24, 0.001e24 );
// MASS_SLIDER.position(width+50, (height-graafikute_kõrgus)*0+320);
// mass_slider_text = createP("Mass:"+mass);
// mass_slider_text.style("font-size","24");
// mass_slider_text.position(width+50, (height-graafikute_kõrgus)*0+270);

// RAADIUS_SLIDER = createSlider(0, 100000, 6371.0, 1);
// RAADIUS_SLIDER.position(width+50,(height-graafikute_kõrgus)*0+400);
// radius_slider_text=createP("Raadius:"+radius);
// radius_slider_text.style("font-size","24");
// radius_slider_text.position(width+50, (height-graafikute_kõrgus)*0+350);

//JUPYTER_button = createButton("Jupyter");
//JUPYTER_button.size(100,50);
//JUPYTER_button.position(width+50, (height-graafikute_kõrgus)*0+420);
//JUPYTER_button.mousePressed(function(){init_model(3);});


// ------------------- PLOT OPTIONS ----------------
  
radio2 = createRadio();
radio2.position(width+50, height/2-90);
radio2.option(1,"Vx(t)");
radio2.option(2,"Vy(t)");
radio2.option(3,"Ax(t)");
radio2.option(4,"Ay(t)");
radio2.option(5,"Etot(t)");
radio2.option(6,"v(|r|)");
radio2.option(7,"a(|r|)");
radio2.style('width', '70px');
radio2.style('checked','0');
radio2.selected('1');

  
PLOT_SCALE_SLIDER=createSlider(0.1, 5, 1, 0.1);
PLOT_SCALE_SLIDER.position(width+50,height-250);
PLOT_SCALE_SLIDER_TEXT=createP("Graafiku skaala");
PLOT_SCALE_SLIDER_TEXT.position(width+50, height-280);
ANIM_SCALE_SLIDER=createSlider(0.01, 5, 1, 0.01);
ANIM_SCALE_SLIDER.position(width+50,height-300);
ANIM_TEXT=createP("Animatsiooni skaala");
ANIM_TEXT.position(width+50,height-330); TIME_STEP_SLIDER=createSlider(0.01,100,10,0.01);
TIME_STEP_SLIDER.position(width+50, height-200);
TIME_STEP_SLIDER_TEXT=createP("Ajasamm h =", h);
TIME_STEP_SLIDER_TEXT.position(width+50,height-230);
//-------------------- PLOT OPTIONS ----------------
  
noLoop();
}


//------------------------------------------ BEGINNING OF DRAW FUNCTION ----------------------------------
function draw() {

clear();
background(0);

  
// mass = MASS_SLIDER.value();
// radius = RAADIUS_SLIDER.value();
plot_scale=PLOT_SCALE_SLIDER.value();
anim_scale=ANIM_SCALE_SLIDER.value();
h=float(TIME_STEP_SLIDER.value());
  
it=it+h;
integration();
  
fill(255);
circle(width/2, (height-graafikute_kõrgus)/2 ,2*radius/(100*anim_scale)); // Planet
  push();
  fill(255,0,0);
circle((yy[0]/(100*anim_scale))+width/2, -1*yy[1]/(100*anim_scale)+(height-graafikute_kõrgus)/2 ,1000/(100*anim_scale)); // Sattelite
  pop();


push();
stroke(0);
strokeWeight(3);
fill(255);
textSize(22);
text("Aeg t: "+ round_3(it),20,30);
text("x:"+round_2(yy[0],2),30,60);
text("y:" + round_2(yy[1],2),30,90);
text("Vx:"+round_2(yy[2],2),30,120);
text("Vy:"+round_2(yy[3],2),30,150);
text("h: "+(round_2(Math.sqrt(yy[0]*yy[0]+yy[1]*yy[1])-radius)), 30,180);
text("Planeet: "+model_name, width/2, (height-graafikute_kõrgus)*0+50);
text(round_0(radius*(anim_scale)),((width-200)+(width-200+radius/100))/2-30,40  );
pop();
// mass_slider_text.html("Mass: "+mass);
// radius_slider_text.html("Radius: "+radius);

  TIME_STEP_SLIDER_TEXT.html("Ajasamm h = "+ h);

// MÕÕTKAVA
push();
stroke(125);
strokeWeight(2);
line(width-200, 50,width-200,75);
line(width-200+radius/100, 50,width-200+radius/100,75);
line(width-200,50,width-200+radius/100,50)
pop();
  
  
  
E_pot = -(G_const*m_body*mass)/(Math.sqrt(yy[0]*yy[0]+yy[1]*yy[1])*1e3);
E_kin = (m_body*(yy[2]*yy[2]+yy[3]*yy[3])*1e6)/2.0;
E_tot = E_pot+E_kin;

r_moodul=Math.sqrt(yy[0]*yy[0]+yy[1]*yy[1]);
v_moodul=Math.sqrt(yy[2]*yy[2]+yy[3]*yy[3]);
a_moodul=Math.sqrt(yp[2]*yp[2]+yp[3]*yp[3]);

// --------------------- PLOT ------------------------------
append(aja_massiiv, it);
append(Vx_massiiv, yy[2]);
append(Vy_massiiv, yy[3]);
append(Ax_massiiv,yp[2]);
append(Ay_massiiv,yp[3]);
append(Epot_massiiv, E_pot);
append(Ekin_massiiv, E_kin);
append(Etot_massiiv, E_tot);  
append(h_massiiv,r_moodul);
append(A_mod_massiiv,a_moodul);
append(V_mod_massiiv,v_moodul);

// plot background
push();
//fill(32,42,68);
fill(25);
rect(0,height-graafikute_kõrgus,width,height);
pop();

// plot xy axis
push();
stroke(255);
strokeWeight(2);
line(plot_offset_X,height-plot_offset_Y,width,height-plot_offset_Y); // X- axis
line(width-10,height-(plot_offset_Y+5),width,height-plot_offset_Y);
line(width-10,height-(plot_offset_Y-5),width,height-plot_offset_Y);
line(plot_offset_X,height,plot_offset_X,height-graafikute_kõrgus); // Y-axis
line(plot_offset_X-5,height-graafikute_kõrgus+10,plot_offset_X, height-graafikute_kõrgus);
line(plot_offset_X+5,height-graafikute_kõrgus+10,plot_offset_X, height-graafikute_kõrgus);
for (j=0; j < 12; j=j+1 ) {
  push();
  stroke(255);
  strokeWeight(1);
  line(plot_offset_X-5,height-25*j, plot_offset_X+5, height-25*j);
  pop();
}
pop();
  
  switch (radio2.value()) {
    //radio value is always a string
    case "1":
      GRAAFIK_Vx_AEG();
      break;
    case "2":
      GRAAFIK_Vy_AEG();
      break;
    case "3":
      GRAAFIK_Ax_AEG();
      break;
    case "4":
      GRAAFIK_Ay_AEG();
      break;
    case "5":
      GRAAFIK_E_tot_AEG();
      break;
    case "6":
      GRAAFIK_Vr();
      break;
    case "7":
      GRAAFIK_Ar();
      break;
    
  }
// ------------------------PLOT-----------------------------
  
if ( Math.sqrt(yy[0]*yy[0]+yy[1]*yy[1])<=radius ) {
  push();
  fill(255,0,0);
  text("Moodul langes planeedi pinnale",width/2,(height-graafikute_kõrgus)/2 );
  pop();
  noLoop();
} 
if (it>nt){
noLoop();
}
}

//------------------------------------------ END OF DRAW FUNCTION ----------------------------------


//------------------------------- RADIO BUTTON PLOTS --------------------------------------------
var V_väärtused=[-15 ,-12.5 ,-10 ,-7.5 ,-5 ,-2.5 ,0 ,2.5 ,5 ,7.5 ,10 ,12.5];

// -------------------------------------------------- Vx(t) plot
function GRAAFIK_Vx_AEG() {
  for (var j=0; j<V_väärtused.length; j=j+1){
    push();
    textAlign(RIGHT);
    text(round_2(V_väärtused[j]*plot_scale) , plot_offset_X-15, height-25*j);
    pop();
  };
  push();
  textAlign(RIGHT);
  textSize(16);
  text("Vx(t)", plot_offset_X-35, height-graafikute_kõrgus/2);
  text("[km/s]", plot_offset_X-35, height-graafikute_kõrgus/2+15);
  text("t [s]",width-25, height-graafikute_kõrgus/2+20 );
  pop();
  push();
  stroke(255,125,0);
  strokeWeight(1);
  for (var k=0; k<=aja_massiiv.length; k=k+1){
    if (abs((Vx_massiiv[k]*10)/plot_scale) <= graafikute_kõrgus/2) {
      point((aja_massiiv[k]/110)+plot_offset_X, -((Vx_massiiv[k]*10)/plot_scale)+(height-plot_offset_Y));
      //console.log(Vx_massiiv[k]);
      if (k>=1) {
        stroke(255,125,0);
        strokeWeight(1);
        line((aja_massiiv[k-1]/110)+plot_offset_X,-Vx_massiiv[k-1]*10/plot_scale+(height-plot_offset_Y) ,(aja_massiiv[k]/110)+plot_offset_X , -Vx_massiiv[k]*10/plot_scale+(height-plot_offset_Y) );
  };
    };
  };
  pop();
};

// -------------------------------------------------- Vy(t) plot
function GRAAFIK_Vy_AEG() {
  for (var j=0; j<V_väärtused.length; j=j+1){
    push();
    textAlign(RIGHT);
    text(round_2(V_väärtused[j]*plot_scale) , plot_offset_X-15, height-25*j);
    pop();
  }
   push();
  textAlign(RIGHT);
  textSize(16);
  text("Vy(t)", plot_offset_X-35, height-graafikute_kõrgus/2);
  text("[km/s]", plot_offset_X-35, height-graafikute_kõrgus/2+15);
  text("t [s]",width-25, height-graafikute_kõrgus/2+20 );
  pop();
push();
stroke(255,125,0);
strokeWeight(1.5);
for (var k=0; k<=aja_massiiv.length; k=k+1) {
      if ( abs((Vy_massiiv[k]*10)/plot_scale) <= graafikute_kõrgus/2 ) {
        point((aja_massiiv[k]/110)+plot_offset_X, -Vy_massiiv[k]*10/plot_scale+(height-plot_offset_Y));
        //console.log(Vy_massiiv[k]);
      if (k>=1) {
        stroke(255,125,0);
        strokeWeight(1);
        line((aja_massiiv[k-1]/110)+plot_offset_X,-Vy_massiiv[k-1]*10/plot_scale+(height-plot_offset_Y) ,(aja_massiiv[k]/110)+plot_offset_X , -Vy_massiiv[k]*10/plot_scale+(height-plot_offset_Y) );
  }
    }
};
pop();
};



// -------------------------------------------------- Ax(t) plot
function GRAAFIK_Ax_AEG(){
  for (var j=0; j<V_väärtused.length; j=j+1){
    push();
    textAlign(RIGHT);
    text(round_4((V_väärtused[j]/100)*plot_scale) , plot_offset_X-15, height-25*j);
    pop();
  };
   push();
  textAlign(RIGHT);
  textSize(16);
  text("Ax(t)", plot_offset_X-35, height-graafikute_kõrgus/2);
  text("[km/s2]", plot_offset_X-35, height-graafikute_kõrgus/2+15);
  text("t [s]",width-25, height-graafikute_kõrgus/2+20 );
  pop();
push();
stroke(255,125,0);
strokeWeight(1.5);
for (var k=0; k<=aja_massiiv.length; k=k+1) {
  if ( abs((Ax_massiiv[k]*1000)/plot_scale) <= graafikute_kõrgus/2 ) {
    point((aja_massiiv[k]/110)+plot_offset_X, -Ax_massiiv[k]*1000/plot_scale+(height-plot_offset_Y));
    //console.log(Ax_massiiv[k]);
      if (k>=1) {
      stroke(255,125,0);
      strokeWeight(1);
      line((aja_massiiv[k-1]/110)+plot_offset_X,-Ax_massiiv[k-1]*1000/plot_scale+(height-plot_offset_Y) ,(aja_massiiv[k]/110)+plot_offset_X , -Ax_massiiv[k]*1000/plot_scale+(height-plot_offset_Y) );
  };
  };
};
pop();    
};

// -------------------------------------------------- Ay(t) plot
function GRAAFIK_Ay_AEG(){
  for (var j=0; j<V_väärtused.length; j=j+1){
    push();
    textAlign(RIGHT);
    text(round_4((V_väärtused[j])/100*plot_scale) , plot_offset_X-15, height-25*j);
    pop();
  };
  push();
  textAlign(RIGHT);
  textSize(16);
  text("Ay(t)", plot_offset_X-35, height-graafikute_kõrgus/2);
  text("[km/s2]", plot_offset_X-35, height-graafikute_kõrgus/2+15);
  text("t [s]",width-25, height-graafikute_kõrgus/2+20 );
  pop();
push();
stroke(255,125,0);
strokeWeight(1.5);
for (var k=0; k<=aja_massiiv.length; k=k+1) {
  if ( abs((Ay_massiiv[k]*1000)/plot_scale) <= graafikute_kõrgus/2 ) {
    point((aja_massiiv[k]/110)+plot_offset_X, -Ay_massiiv[k]*1000/plot_scale+(height-plot_offset_Y));
    //console.log(Ay_massiiv[k]);
    if (k>=1) {
      stroke(255,125,0);
      strokeWeight(1);
      line((aja_massiiv[k-1]/110)+plot_offset_X,-Ay_massiiv[k-1]*1000/plot_scale+(height-plot_offset_Y) ,(aja_massiiv[k]/110)+plot_offset_X , -Ay_massiiv[k]*1000/plot_scale+(height-plot_offset_Y) );
  };
  };
};
pop();    
};

//-------------------------------------------------- Etot(t) plot
function GRAAFIK_E_tot_AEG(){
    for (var j=0; j<V_väärtused.length; j=j+1){
    push();
    textAlign(RIGHT);
    text(round_4((V_väärtused[j])*10*plot_scale) , plot_offset_X-15, height-25*j);
    pop();
  };
  push();
  textAlign(RIGHT);
  textSize(16);
  text("E(t)", plot_offset_X-35, height-graafikute_kõrgus/2);
  text("[J]", plot_offset_X-35, height-graafikute_kõrgus/2+15);
  text("t [s]",width-25, height-graafikute_kõrgus/2+20 );
  pop();
  push();
  textAlign(LEFT);
  textSize(14);
  text("1e8", plot_offset_X+10, height-graafikute_kõrgus+15);
  pop();
  
push();
stroke(255,125,0);
strokeWeight(1.5);
for (var k=0; k<=aja_massiiv.length; k=k+1) {
  if (abs((Etot_massiiv[k]/1e8)/plot_scale) <= graafikute_kõrgus/2 ) {
    push();
    fill(255);
    stroke(255,125,0);
    strokeWeight(1);
    point((aja_massiiv[k]/110)+plot_offset_X, -(Etot_massiiv[k]/1e8)/plot_scale+(height-plot_offset_Y));
    //console.log(Etot_massiiv[k]/1e8);
    pop();
    if (k>=1) {
      stroke(255,125,0);
      strokeWeight(1);
      line((aja_massiiv[k-1]/110)+plot_offset_X, -(Etot_massiiv[k-1]/1e8)/plot_scale+(height-plot_offset_Y) ,(aja_massiiv[k]/110)+plot_offset_X , -(Etot_massiiv[k]/1e8)/plot_scale+(height-plot_offset_Y) );
  };
  };
};
pop();      
};


// -------------------------------------------------- V(|r|) plot
function GRAAFIK_Vr(){
    for (var j=0; j<V_väärtused.length; j=j+1){
      push();
      textAlign(RIGHT);
      text(round_4((V_väärtused[j])*10*plot_scale) , plot_offset_X-15, height-25*j);
      pop();
  };
  push();
  textAlign(RIGHT);
  textSize(16);
  text("V(|r|)", plot_offset_X-35, height-graafikute_kõrgus/2);
  text("[km/s]", plot_offset_X-35, height-graafikute_kõrgus/2+15);
  text("|r| [km]",width-25, height-graafikute_kõrgus/2+20 );
  pop();
push();
stroke(255,125,0);
strokeWeight(1.5);
for (var k=0; k<=h_massiiv.length; k=k+1) {
  if (abs((V_mod_massiiv[k])/plot_scale) <= graafikute_kõrgus/2 ) {
    push();
    fill(255);
    stroke(255,125,0);
    strokeWeight(1);
    point((h_massiiv[k]/110)+plot_offset_X, -(V_mod_massiiv[k])/plot_scale+(height-plot_offset_Y));
    //console.log(V_mod_massiiv[k]);
    pop();
    if (k>=1) {
      stroke(255,125,0);
      strokeWeight(1);
      line((h_massiiv[k-1]/110)+plot_offset_X, -(V_mod_massiiv[k-1])/plot_scale+(height-plot_offset_Y) ,(h_massiiv[k]/110)+plot_offset_X , -(V_mod_massiiv[k])/plot_scale+(height-plot_offset_Y) );
  };
  };
};
pop();      
};


// -------------------------------------------------- A(|r|) plot
function GRAAFIK_Ar(){
      for (var j=0; j<V_väärtused.length; j=j+1){
      push();
      textAlign(RIGHT);
      text(round_4((V_väärtused[j])/100*plot_scale) , plot_offset_X-15, height-25*j);
      pop();
  };
  push();
  textAlign(RIGHT);
  textSize(16);
  text("A(|r|)", plot_offset_X-35, height-graafikute_kõrgus/2);
  text("[km/s2]", plot_offset_X-35, height-graafikute_kõrgus/2+15);
  text("|r| [km]",width-25, height-graafikute_kõrgus/2+20 );
  pop();
push();
stroke(255,125,0);
strokeWeight(1.5);
for (var k=0; k<=h_massiiv.length; k=k+1) {
  if (abs((A_mod_massiiv[k]*1e3)/plot_scale) <= graafikute_kõrgus/2 ) {
    push();
    fill(255);
    stroke(255,125,0);
    strokeWeight(1);
    point((h_massiiv[k]/110)+plot_offset_X, -(A_mod_massiiv[k]*1e3)/plot_scale+(height-plot_offset_Y));
    //console.log(A_mod_massiiv[k]);
    pop();
    if (k>=1) {
      stroke(255,125,0);
      strokeWeight(1);
      line((h_massiiv[k-1]/110)+plot_offset_X, -(A_mod_massiiv[k-1]*1e3)/plot_scale+(height-plot_offset_Y) ,(h_massiiv[k]/110)+plot_offset_X , -(A_mod_massiiv[k]*1e3)/plot_scale+(height-plot_offset_Y) );
  };
  };
};
pop(); 
};



// ------------------------------- PROPER ARITHMETIC ROUNDING FUNCTIONS --------------------------------
function round_0(v) {
  return (Math.sign(v)*Math.round(Math.abs(v)));
}

function round_1(v) {
  return (Math.sign(v)*Math.round(Math.abs(v)*10)/10);
}

function round_2(v) {
  return (Math.sign(v)*Math.round(Math.abs(v)*100)/100);
}

function round_3(v) {
  return (Math.sign(v)*Math.round(Math.abs(v)*1000)/1000);
}

function round_4(v) {
  return (Math.sign(v)*Math.round(Math.abs(v)*10000)/10000);
}
