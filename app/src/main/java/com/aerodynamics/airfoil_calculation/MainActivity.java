package com.aerodynamics.airfoil_calculation;

////////////////////////////////////////////////////////////////////////////////////////////////////
//Questa classe contiene tutto il codice Java per la schermata iniziale dell'applicazione	      //
////////////////////////////////////////////////////////////////////////////////////////////////////

import androidx.appcompat.app.AppCompatActivity;
import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.EditText;
import android.widget.Toast;
import com.google.firebase.analytics.FirebaseAnalytics;

public class MainActivity extends AppCompatActivity {

	//Definisco alcune variabili che verranno inizializzate in seguito.
	//Tipo EditText --> E' il tipo di variabile corrispondente al campo di testo modificabile dall'utente
	//Le variabili di tipo String, double, int sono quelle in cui verranno salvati i valori immessi
	//dall'utente.

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

	//le seguenti 2 righe vengono inserite di default per ogni activity, il loro significato è a me
	//abbastanza oscuro anche se è possibile trovare facilmente una spiegazione su internet.
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);

		//Indico quale file di layout .xml corrisponde a questa activity
		setContentView(R.layout.activity_main);

		//Firebase è un prodotto di Google che ha molte funzionalità, tra le quali l'analisi
		//delle statistiche dell'utilizzo dell'app da parte degli utenti.
		mFirebaseAnalytics = FirebaseAnalytics.getInstance(this);

		//findViewById permette di associare una variabile (ad esempio di tipo EditText) con la
		//corrispondente porzione di codice all'interno del file di layout .xml
		nacaDes = findViewById(R.id.addNACAd);
		//setHint permette di fornire un suggerimento per l'immisione del testo
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

	//Nel file activity_main ho aggiunto agli oggetti di tipo "button" la
	//riga di codice android:onClick="onClick". Ciò significa che nel caso in cui l'utente ci clicchi
	//sopra, questi oggetti sono predisposti a fare accadere qualcosa. Di seguito ho scritto il metodo
	//onClick che definisce cosa succede quando l'utente clicca su un oggetto di tipo "Button"


	public void onClick(View view) {

		//la funzione getText permette di ottenere il testo immesso dall'utente. toString è necessario
		//per convertire il tipo EditText in String
		nacaD = nacaDes.getText().toString();

		//Per il debug, senza stare a scrivere qualcosa ogni volta ho messo dei controlli per ogni
		//testo da inserire. Nel caso il campo di testo sia vuoto (isEmpty), assegno di default
		//dei valori alle variabili corrispondenti, per risparmiare il tempo di immissione dei caratteri.

		if(nacaD.isEmpty()){
			nacaD = "0012";
		}
		if(chordLen.getText().toString().isEmpty()){
			chLn = 1.0;
		}
		else

		//La funzione parseDouble permette di convertire la stringa di testo recuperata dalla EditText
		//attraverso getText().toString() in un valore di tipo double.
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

		//Se il profilo scelto è a 4 o 5 cifre, e contemporaneamente il numero di pannelli
		//selezionato è compreso tra 10 e 200, allora tutto è ok
		if((nacaD.length() == 4 || nacaD.length() == 5) && (nPan > 9 && nPan < 201)) {

			//Definisco nello specifico cosa fare nel caso in cui venga selezionato l'oggetto del
			//layout button1
			if(view.getId() == R.id.button1) {

				//Con le righe di codice seguenti definisco che se viene selezionato il primo button,
				//l'activity che comparirà a schermo sarà quella corrispondente alla classe BuildGeometry

				//Attraverso il metodo putExtra posso chiedere che da un'activity a quell'altra vengano
				//trasmessi i valori desiderati di alcune variabili, nello specifico quelle corrispondenti
				//alle scelte dell'utente.
				Intent sendIntent = new Intent(MainActivity.this, BuildGeometry.class);
				sendIntent.putExtra("NACA_D", nacaD);
				sendIntent.putExtra("Chord", chLn);
				sendIntent.putExtra("Stagger", incid);
				sendIntent.putExtra("nPan", nPan);
				sendIntent.putExtra("Pressure", pr);
				sendIntent.putExtra("Density", den);
				sendIntent.putExtra("Speed", spd);
				sendIntent.putExtra("AngleOA", aOA);

				//startActivity fa esattamente quello che dice il nome. In questo caso fa partire
				//l'activity corrispondente a quanto definito dall'oggetto di tipo Intent sendIntent,
				//ossia l'activity BuildGeometry.
				startActivity(sendIntent);
			}
			else{

				//Caso button2
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

		// Di seguito definisco cosa fare se quanto immesso dall'utente non soddisfi i vincoli presenti
		else if (nacaD.length() < 4 || nacaD.length() > 5){

			//Toast è un oggetto che permette di stampare a schermo un messaggio per la durata di qualche secondo,
			//in questo caso spiega all'utente cosa ha sbagliato e perchè non passa alla schermata successiva.

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


