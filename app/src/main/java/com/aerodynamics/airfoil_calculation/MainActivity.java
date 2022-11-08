package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
// This class contains all the Java code about the first and main screen of the app      	      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import androidx.appcompat.app.AppCompatActivity;
import androidx.constraintlayout.widget.ConstraintLayout;

import android.content.Intent;
import android.graphics.Color;
import android.graphics.PorterDuff;
import android.os.Bundle;
import android.view.ContextMenu;
import android.view.Menu;
import android.view.MenuInflater;
import android.view.MenuItem;
import android.view.View;
import android.widget.Button;
import android.widget.CompoundButton;
import android.widget.EditText;
import android.widget.Toast;

import com.google.android.material.switchmaterial.SwitchMaterial;
import com.google.firebase.analytics.FirebaseAnalytics;

public class MainActivity extends AppCompatActivity {

	// In the following lines of code some variable are defined, which later will be initialized.
	// The variable type SwitchMaterial is the widget which correspond to the "True/False" choice.
	// The variable type EditText is the widget which correspond to the text field which can be filled by the user.
	// The variable types String, double, int are the ones in which will be saved the numbers written
	// by the user.

	SwitchMaterial groundEffect, wingTail;
	EditText nacaDesW, chordLenW, incidenceW, nPanels, speed, anglofattack, groundHeight, nacaDesT,
			 chordLenT, incidenceT, wingTailX, wingTailY;
	String nacaDW, nacaDT;
	double chLnW, incidW, spd, aOA, grndH, chLnT, incidT, wTX, wTY;
	int nPan;
	boolean isGround, isWT;
	//EditText pressure, density;
	//double pr, den;

	@Override

	// The following two lines of code are written by default in every activity. Their meaning is not
	// so clear to me, anyway it is explained on the suggested book.
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);

		// This line sets the correspondent layout file (.xml) to this activity
		setContentView(R.layout.activity_main);

		// Firebase is a Google service which has many features, among them there is the possibility
		// to analise many statistics about the user behaviour inside th app.
		// With the following line of code the service will start to collect the default statistics.
		FirebaseAnalytics mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		// findViewById allows to associate one variable (for example an EditText type one)
		// with the correspondent portion of code inside the layout (.xml) file.

		Button button = findViewById(R.id.button11);
		int defCol = button.getTextColors().getDefaultColor();
		Button button1 = findViewById(R.id.button3);
		Button button4 = findViewById(R.id.button4); //Generate dataset button
		Button button5 = findViewById(R.id.button5);

		groundEffect = findViewById(R.id.switchGround);
		groundEffect.setOnCheckedChangeListener((buttonView, isChecked) -> {
			// If the user select the ground effect option, the layout relative to the ground
			// effect details will be displayed.
			if (isChecked) {
				ConstraintLayout constraintLayout = findViewById(R.id.all05);
				constraintLayout.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout1 = findViewById(R.id.all5);
				constraintLayout1.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout2 = findViewById(R.id.s9);
				constraintLayout2.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout3 = findViewById(R.id.d9);
				constraintLayout3.setVisibility(View.VISIBLE);
				isGround = true;
				button.setEnabled(false);
				button.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
				button.setTextColor(Color.LTGRAY);
				button4.setEnabled(false);
				button4.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
				button4.setTextColor(Color.LTGRAY);
				button5.setEnabled(false);
				button5.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
				button5.setTextColor(Color.LTGRAY);
			} else {
				ConstraintLayout constraintLayout = findViewById(R.id.all05);
				constraintLayout.setVisibility(View.GONE);
				ConstraintLayout constraintLayout1 = findViewById(R.id.all5);
				constraintLayout1.setVisibility(View.GONE);
				ConstraintLayout constraintLayout2 = findViewById(R.id.s9);
				constraintLayout2.setVisibility(View.GONE);
				ConstraintLayout constraintLayout3 = findViewById(R.id.d9);
				constraintLayout3.setVisibility(View.GONE);
				button.setEnabled(true);
				button.getBackground().setColorFilter(null);
				button.setTextColor(defCol);
				button4.setEnabled(true);
				button4.getBackground().setColorFilter(null);
				button4.setTextColor(defCol);
				button5.setEnabled(true);
				button5.getBackground().setColorFilter(null);
				button5.setTextColor(defCol);
				isGround = false;
			}
		});
		wingTail = findViewById(R.id.switchTail);
		wingTail.setOnCheckedChangeListener((buttonView, isChecked) -> {
			// If the user select the ground effect option, the layout relative to the ground
			// effect details will be displayed.
			if (isChecked) {
				ConstraintLayout constraintLayout = findViewById(R.id.all06);
				constraintLayout.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout1 = findViewById(R.id.all6);
				constraintLayout1.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout2 = findViewById(R.id.s10);
				constraintLayout2.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout3 = findViewById(R.id.d10);
				constraintLayout3.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout4 = findViewById(R.id.s11);
				constraintLayout4.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout5 = findViewById(R.id.d11);
				constraintLayout5.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout6 = findViewById(R.id.s12);
				constraintLayout6.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout7 = findViewById(R.id.d12);
				constraintLayout7.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout8 = findViewById(R.id.s13);
				constraintLayout8.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout9 = findViewById(R.id.d13);
				constraintLayout9.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout10 = findViewById(R.id.s14);
				constraintLayout10.setVisibility(View.VISIBLE);
				ConstraintLayout constraintLayout11 = findViewById(R.id.d14);
				constraintLayout11.setVisibility(View.VISIBLE);
				button.setEnabled(false);
				button1.setEnabled(false);
				button.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
				button1.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
				button.setTextColor(Color.LTGRAY);
				button1.setTextColor(Color.LTGRAY);
				button4.setEnabled(false);
				button4.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
				button4.setTextColor(Color.LTGRAY);
				button5.setEnabled(false);
				button5.getBackground().setColorFilter(Color.rgb(245,245,245), PorterDuff.Mode.MULTIPLY);
				button5.setTextColor(Color.LTGRAY);
				isWT = true;
			} else {
				ConstraintLayout constraintLayout = findViewById(R.id.all06);
				constraintLayout.setVisibility(View.GONE);
				ConstraintLayout constraintLayout1 = findViewById(R.id.all6);
				constraintLayout1.setVisibility(View.GONE);
				ConstraintLayout constraintLayout2 = findViewById(R.id.s10);
				constraintLayout2.setVisibility(View.GONE);
				ConstraintLayout constraintLayout3 = findViewById(R.id.d10);
				constraintLayout3.setVisibility(View.GONE);
				ConstraintLayout constraintLayout4 = findViewById(R.id.s11);
				constraintLayout4.setVisibility(View.GONE);
				ConstraintLayout constraintLayout5 = findViewById(R.id.d11);
				constraintLayout5.setVisibility(View.GONE);
				ConstraintLayout constraintLayout6 = findViewById(R.id.s12);
				constraintLayout6.setVisibility(View.GONE);
				ConstraintLayout constraintLayout7 = findViewById(R.id.d12);
				constraintLayout7.setVisibility(View.GONE);
				ConstraintLayout constraintLayout8 = findViewById(R.id.s13);
				constraintLayout8.setVisibility(View.GONE);
				ConstraintLayout constraintLayout9 = findViewById(R.id.d13);
				constraintLayout9.setVisibility(View.GONE);
				ConstraintLayout constraintLayout10 = findViewById(R.id.s14);
				constraintLayout10.setVisibility(View.GONE);
				ConstraintLayout constraintLayout11 = findViewById(R.id.d14);
				constraintLayout11.setVisibility(View.GONE);
				button.setEnabled(true);
				button1.setEnabled(true);
				button.getBackground().setColorFilter(null);
				button1.getBackground().setColorFilter(null);
				button.setTextColor(defCol);
				button1.setTextColor(defCol);
				button4.setEnabled(true);
				button4.getBackground().setColorFilter(null);
				button4.setTextColor(defCol);
				button5.setEnabled(true);
				button5.getBackground().setColorFilter(null);
				button5.setTextColor(defCol);
				isWT = false;
			}
		});

		nacaDesW = findViewById(R.id.addnacaDW);
		// setHint allows to display an hint for what the user can type inside the text field.
		nacaDesW.setHint("0012");
		chordLenW = findViewById(R.id.addChordW);
		chordLenW.setHint("1.0");
		incidenceW = findViewById(R.id.editStaggerW);
		incidenceW.setHint("0.0");
		nPanels = findViewById(R.id.editNPan);
		nPanels.setHint("20");
		speed = findViewById(R.id.editSpeed);
		speed.setHint("30.0");
		anglofattack = findViewById(R.id.editAngOAtt);
		anglofattack.setHint("0.0");
		//pressure = findViewById(R.id.addPressure);
		//pressure.setHint("101325");
		//density = findViewById(R.id.addDensity);
		//density.setHint("1.225");
		groundHeight = findViewById(R.id.editGroundHeight);
		groundHeight.setHint("1.0");
		nacaDesT = findViewById(R.id.addnacaDT);
		nacaDesT.setHint("0012");
		chordLenT = findViewById(R.id.addChordT);
		chordLenT.setHint("1.0");
		incidenceT = findViewById(R.id.editStaggerT);
		incidenceT.setHint("0.0");
		wingTailX = findViewById(R.id.editWingTailX);
		wingTailX.setHint("1.5");
		wingTailY = findViewById(R.id.editWingTailY);
		wingTailY.setHint("0.2");
	}

	// Inside the layout file activity_main some "button" widgets type has been defined with the
	// following attribute: android:onClick="onClick". This attribute makes literally the button clickable,
	// which means that these widgets are allowed to make something happen once they are clicked.
	// On the following lines of code the onClick method has been defined, which contains the instruction
	// that needs to be executed once the buttons are clicked.

	public void onClick(View view) {
		// The getText method returns the text filled by the user. The toString method is necessary
		// to convert the EditText type to a string.
		nacaDW = nacaDesW.getText().toString();
		// To speed up the debug, without the need to type some values every time, I've added
		// the following if statements which will assign a value to the variables if
		// the corresponding text fields are empty.
		if(nacaDW.isEmpty()){
			nacaDW = "0012";
		}

		if(chordLenW.getText().toString().isEmpty()){
			chLnW = 1.0;
		}
		else
		// The parseDouble method allows to convert the text string retrieved by the getText().toString()
		// methods into a double type value.
		chLnW = Double.parseDouble(chordLenW.getText().toString());

		if(incidenceW.getText().toString().isEmpty()){
			incidW = 0.0;
		}
		else
		incidW = Double.parseDouble(incidenceW.getText().toString());

		if(nPanels.getText().toString().isEmpty()){
			nPan = 20;
		}
		else
		nPan = Integer.parseInt(nPanels.getText().toString());

		//if(pressure.getText().toString().isEmpty()){
		//	pr = 101325.0;
		//}
		//else
		//pr = Double.parseDouble(pressure.getText().toString());
//
		//if(density.getText().toString().isEmpty()){
		//	den = 1.225;
		//}
		//else
		//den = Double.parseDouble(density.getText().toString());

		if(speed.getText().toString().isEmpty()){
			spd = 30.0;
		}
		else
		spd = Double.parseDouble(speed.getText().toString());

		if(anglofattack.getText().toString().isEmpty()){
			aOA = 0.0;
		}
		else
		aOA = Double.parseDouble(anglofattack.getText().toString());

		if(groundHeight.getText().toString().isEmpty()){
			grndH = 1.0;
		}
		else
		grndH = Double.parseDouble(groundHeight.getText().toString());

		nacaDT = nacaDesT.getText().toString();

		if(nacaDT.isEmpty()){
			nacaDT = "0012";
		}

		if(chordLenT.getText().toString().isEmpty()){
			chLnT = 1.0;
		}
		else
			chLnT = Double.parseDouble(chordLenT.getText().toString());

		if(incidenceT.getText().toString().isEmpty()){
			incidT = 0.0;
		}
		else
			incidT = Double.parseDouble(incidenceT.getText().toString());

		if(wingTailX.getText().toString().isEmpty()){
			wTX = 1.5;
		}
		else
		wTX = Double.parseDouble((wingTailX.getText().toString()));

		if(wingTailY.getText().toString().isEmpty()){
			wTY = 0.2;
		}
		else
		wTY = Double.parseDouble((wingTailY.getText().toString()));

		// If the airfoil chosen is characterized by 4 or 5 digits, and at the same time the number
		// of panels selected is between 10 and 200, everything is fine. For 5 digits airfoil the
		// designation, in addition, can start only with number 210, 220, 230, 240, 250.
		if((nacaDW.length() == 4 || (((nacaDW.length() == 5)) && ((nacaDW.startsWith("210")||
				nacaDW.startsWith("220")||nacaDW.startsWith("230")
				||nacaDW.startsWith("240")||nacaDW.startsWith("250"))))) && (nPan > 9 && nPan < 201)
		&&((nacaDT.length() == 4 || (((nacaDT.length() == 5)) && ((nacaDT.startsWith("210")||
				nacaDT.startsWith("220")||nacaDT.startsWith("230")
				||nacaDT.startsWith("240")||nacaDT.startsWith("250"))))))) {

			// Now it's defined what happens when the widget button1 is clicked
			if(view.getId() == R.id.button1) {

				// Through an Intent, once the button1 is clicked, the BuildGeometry activity will
				// be displayed on the device sceen.

				// Thanks to the putExtra method some data can be transferred between the two activities,
				// in the specific the ones corresponding to the data filled by the user.
				Intent sendIntent = new Intent(MainActivity.this, BuildGeometry.class);
				sendIntent.putExtra("NACA_DW", nacaDW);
				sendIntent.putExtra("ChordW", chLnW);
				sendIntent.putExtra("StaggerW", incidW);
				sendIntent.putExtra("nPan", nPan);
				sendIntent.putExtra("GroundH", grndH);
				sendIntent.putExtra("NACA_DT", nacaDT);
				sendIntent.putExtra("ChordT", chLnT);
				sendIntent.putExtra("StaggerT", incidT);
				sendIntent.putExtra("checkGround", isGround);
				sendIntent.putExtra("checkTail", isWT);
				sendIntent.putExtra("wTX", wTX);
				sendIntent.putExtra("wTY", wTY);

				// The startActivity method does exactly what the name says.
				startActivity(sendIntent);
			}
			// button2 case
			else if (view.getId() == R.id.button2) {
				Intent sendIntent = new Intent(MainActivity.this, BuildCoefficient.class);
				sendIntent.putExtra("NACA_DW", nacaDW);
				sendIntent.putExtra("ChordW", chLnW);
				sendIntent.putExtra("StaggerW", incidW);
				sendIntent.putExtra("nPan", nPan);
				//sendIntent.putExtra("Pressure", pr);
				//sendIntent.putExtra("Density", den);
				sendIntent.putExtra("Speed", spd);
				sendIntent.putExtra("AngleOA", aOA);
				sendIntent.putExtra("GroundH", grndH);
				sendIntent.putExtra("NACA_DT", nacaDT);
				sendIntent.putExtra("ChordT", chLnT);
				sendIntent.putExtra("StaggerT", incidT);
				sendIntent.putExtra("WingTailX", wTX);
				sendIntent.putExtra("WingTailY", wTY);
				sendIntent.putExtra("checkGround", isGround);
				sendIntent.putExtra("checkTail", isWT);

				startActivity(sendIntent);
			}
			// button11 case
			else if(view.getId() == R.id.button11){
				Intent sendIntent = new Intent(MainActivity.this, ThinAirfoil.class);
				sendIntent.putExtra("NACA_DW", nacaDW);
				sendIntent.putExtra("ChordW", chLnW);
				sendIntent.putExtra("StaggerW", incidW);
				sendIntent.putExtra("nPan", nPan);
				sendIntent.putExtra("AngleOA", aOA);

				startActivity(sendIntent);
			}
			else if(view.getId() == R.id.button3){

				// button3 case
				Intent sendIntent = new Intent(MainActivity.this, ClAlpha.class);
				sendIntent.putExtra("NACA_DW", nacaDW);
				sendIntent.putExtra("ChordW", chLnW);
				sendIntent.putExtra("nPan", nPan);
				//sendIntent.putExtra("Pressure", pr);
				//sendIntent.putExtra("Density", den);
				sendIntent.putExtra("Speed", spd);
				sendIntent.putExtra("GroundH", grndH);
				sendIntent.putExtra("checkGround", isGround);

				startActivity(sendIntent);
			}
			else if(view.getId() == R.id.button4){
				Intent sendIntent = new Intent(MainActivity.this, GenerateDataset.class);
				startActivity(sendIntent);
			}
			else{
				Intent sendIntent = new Intent(MainActivity.this, AI_Inverse_problem.class);
				startActivity(sendIntent);
			}
		}

		// In the following lines it's defined what happens when the user types something which
		// does not satisfy the constrains imposed by me
		else if (nacaDW.length() < 4 || nacaDW.length() > 5 || nacaDT.length() < 4 || nacaDT.length() > 5){

			// The Toast is an object which prints on the screen a message for few seconds,
			// in this case it explains to the user what went wrong.

			Toast.makeText(this, R.string.naca45Error, Toast.LENGTH_LONG).show();
		}
		else if (nPan <= 9){
			Toast.makeText(this, R.string.nPanErrorMin, Toast.LENGTH_LONG).show();
		}
		else if ((nacaDW.length() == 5 && !(nacaDW.startsWith("210")||
				nacaDW.startsWith("220")||nacaDW.startsWith("230")
				||nacaDW.startsWith("240")||nacaDW.startsWith("250")))||
				(nacaDT.length() == 5 && !(nacaDT.startsWith("210")||
						nacaDT.startsWith("220")||nacaDT.startsWith("230")
						||nacaDT.startsWith("240")||nacaDT.startsWith("250")))){
			Toast.makeText(this, R.string.naca5Error, Toast.LENGTH_LONG).show();
		}
		else {
			Toast.makeText(this, R.string.nPanErrorMax, Toast.LENGTH_LONG).show();
		}
	}
//	@Override
//	public boolean onCreateOptionsMenu(Menu menu) {
//		MenuInflater inflater = getMenuInflater();
//		inflater.inflate(R.menu.my_menu, menu);
//		return true;
//	}
//	@Override
//	public boolean onOptionsItemSelected(MenuItem item) {
//		if (item.getItemId() == R.id.info) {
//			Intent intent = new Intent(MainActivity.this, InfoActivity.class);
//			startActivity(intent);
//			return true;
//		}
//		return super.onContextItemSelected(item);
//	}
}


