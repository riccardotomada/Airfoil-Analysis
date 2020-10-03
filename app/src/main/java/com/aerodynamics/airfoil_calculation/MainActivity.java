package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
// This class contains all the Java code about the first and main screen of the app      	      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import androidx.appcompat.app.AppCompatActivity;
import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.EditText;
import android.widget.Toast;
import com.google.firebase.analytics.FirebaseAnalytics;

public class MainActivity extends AppCompatActivity {

	// In the following lines of code some variable are defined, which later will be initialized.
	// The variable type EditText is the widget which correspond to the text field which can be filled by the user.
	// The variable types String, double, int are the ones in which will be saved the numbers written
	// by the user.

	EditText nacaDes;
	String nacaD;
	EditText chordLen;
	double chLn;
	EditText incidence;
	double incid;
	EditText nPanels;
	int nPan;
	EditText pressure;
	double pr;
	EditText density;
	double den;
	EditText speed;
	double spd;
	EditText anglofattack;
	double aOA;

	private FirebaseAnalytics mFirebaseAnalytics;

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
		mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		// findViewById allows to associate one variable (for example an EditText type one)
		// with the correspondent portion of code inside the layout (.xml) file.
		nacaDes = findViewById(R.id.addNACAd);
		// setHint allows to display an hint for what the user can type inside the text field.
		nacaDes.setHint("0012");
		chordLen = findViewById(R.id.addChord);
		chordLen.setHint("1.0");
		incidence = findViewById(R.id.editStagger);
		incidence.setHint("0.0");
		nPanels = findViewById(R.id.editNPan);
		nPanels.setHint("20");
		speed = findViewById(R.id.editSpeed);
		speed.setHint("30.0");
		anglofattack = findViewById(R.id.editAngOAtt);
		anglofattack.setHint("0.0");
		pressure = findViewById(R.id.addPressure);
		pressure.setHint("101325");
		density = findViewById(R.id.addDensity);
		density.setHint("1.225");
	}

	// Inside the layout file activity_main some "button" widgets type has been defined with the
	// following attribute: android:onClick="onClick". This attribute makes literally the button clickable,
	// which means that these widgets are allowed to make something happen once they are clicked.
	// On the following lines of code the onClick method has been defined, which contains the instruction
	// that needs to be executed once the buttons are clicked.


	public void onClick(View view) {

		// The getText method returns the text filled by the user. The toString method is necessary
		// to convert the EditText type to a string.
		nacaD = nacaDes.getText().toString();

		// To speed up the debug, without the need to type some values every time, I've added
		// the following if statements which will assign a value to the variables if
		// the corresponding text fields are empty.

		if(nacaD.isEmpty()){
			nacaD = "0012";
		}
		if(chordLen.getText().toString().isEmpty()){
			chLn = 1.0;
		}
		else

		// The parseDouble method allows to convert the text string retrieved by the getText().toString()
		// methods into a double type value.
		chLn = Double.parseDouble(chordLen.getText().toString());

		if(incidence.getText().toString().isEmpty()){
			incid = 0.0;
		}
		else
		incid = Double.parseDouble(incidence.getText().toString());

		if(nPanels.getText().toString().isEmpty()){
			nPan = 20;
		}
		else
		nPan = Integer.parseInt(nPanels.getText().toString());

		if(pressure.getText().toString().isEmpty()){
			pr = 101325.0;
		}
		else
		pr = Double.parseDouble(pressure.getText().toString());

		if(density.getText().toString().isEmpty()){
			den = 1.225;
		}
		else
		den = Double.parseDouble(density.getText().toString());

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

		// If the airfoil chosen is characterized by 4 or 5 digits, and at the same time the number
		// of panels selected is between 10 and 200, everything is fine
		if((nacaD.length() == 4 || nacaD.length() == 5) && (nPan > 9 && nPan < 201)) {

			// Now it's defined what happens when the widget button1 is clicked
			if(view.getId() == R.id.button1) {

				// Through an Intent, once the button1 is clicked, the BuildGeometry activity will
				// be displayed on the device sceen.

				// Thanks to the putExtra method some data can be transferred between the two activities,
				// in the specific the ones corresponding to the data filled by the user.
				Intent sendIntent = new Intent(MainActivity.this, BuildGeometry.class);
				sendIntent.putExtra("NACA_D", nacaD);
				sendIntent.putExtra("Chord", chLn);
				sendIntent.putExtra("Stagger", incid);
				sendIntent.putExtra("nPan", nPan);
				sendIntent.putExtra("Pressure", pr);
				sendIntent.putExtra("Density", den);
				sendIntent.putExtra("Speed", spd);
				sendIntent.putExtra("AngleOA", aOA);

				// The startActivity method does exactly what the name says.
				startActivity(sendIntent);
			}
			else{

				// button2 case
				Intent sendIntent = new Intent(MainActivity.this, ClAlpha.class);
				sendIntent.putExtra("NACA_D", nacaD);
				sendIntent.putExtra("Chord", chLn);
				sendIntent.putExtra("nPan", nPan);
				sendIntent.putExtra("Pressure", pr);
				sendIntent.putExtra("Density", den);
				sendIntent.putExtra("Speed", spd);
				startActivity(sendIntent);
			}
		}

		// In the following lines it's defined what happens when the user types something which
		// does not satisfy the constrains imposed by me
		else if (nacaD.length() < 4 || nacaD.length() > 5){

			// The Toast is an object which prints on the screen a message for few seconds,
			// in this case it explains to the user what went wrong.

			Toast.makeText(this, "Only NACA 4 and 5 digit available", Toast.LENGTH_LONG).show();
		}
		else if (nPan <= 9){
			Toast.makeText(this, "Too few Panels chosen for the discretization", Toast.LENGTH_LONG).show();
		}
		else {
			Toast.makeText(this, "Too many Panels chosen for the discretization", Toast.LENGTH_LONG).show();
		}

	}
}


