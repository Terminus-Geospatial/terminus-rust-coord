
use std::env;
use std::collections::LinkedList;


pub struct CmdArgs {

}

fn get_value( x: Option<String> ) -> String {

    if let Some( val ) = x {
        return val;
    }
    panic!("Failed to parse next field.");
}

/**
 * Simple Command-Line application for performing coordinate conversions.
 */
 pub fn parse_command_line() -> CmdArgs {

    // Collect arguments to be parsed
    let mut args: LinkedList<String> = env::args().collect();

    // Copies of arguments
    let mut input_coord: Vec<f64> = Vec::new();

    // Fetch all parameters
    while args.len() > 0 {

        // Fetch the string out of the linked list
        let arg = args.pop_front();

        if let Some(a) = arg {
            
            let flag = String::from(a);
            
            // Input Coordinate CSV path
            if flag == "-i" {
                
                while args.len() > 0 && get_value(args.front().cloned()).starts_with('-') {

                    // Fetch the next argument 
                    let next_val: f64 = get_value( args.pop_front() ).parse().unwrap();
                    input_coord.push( next_val );
                }
                
            }

        } else {
            println!("error fetching parameter");
        }

    }

    // Construct return object
    let config = CmdArgs{
        
    };

    return config;
}